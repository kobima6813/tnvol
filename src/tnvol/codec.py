# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Yackov Mozgori
"""tnvol codec: 3D CDF 9/7 wavelet + packet selection + dead-zone RDO +
context-mixing range coder + Wiener post-filter + L-infinity correction +
palette mode.

Public helpers
--------------
compress(volume, *, psnr=None, linf=None, threads=1) -> bytes
decompress(blob) -> (array, tau)
info(blob) -> dict with shape, dtype, tau

The native library is loaded from the built ``tnvol._libwlc`` extension, or
from a ``libwlc.so`` sitting next to this file.
"""
from __future__ import annotations

import ctypes
import glob
import os
import struct
import sys
import tempfile

# Defaults match the frozen r10 encoder. Override with WLC_* environment variables.
KT = float(os.environ.get("WLC_KT", 0.0))
RND = float(os.environ.get("WLC_RND", 0.38))
MAXL = int(os.environ.get("WLC_MAXL", 6))
LAM = float(os.environ.get("WLC_LAM", 0.08))
KTL = float(os.environ.get("WLC_KTL", -1))
FK = int(os.environ.get("WLC_FK", 6))

P = ctypes.c_void_p
I = ctypes.c_int
Lg = ctypes.c_long
D = ctypes.c_double


def _load_lib():
    candidates = []
    try:
        from importlib.util import find_spec
        spec = find_spec("tnvol._libwlc")
        if spec is not None and getattr(spec, "origin", None) and os.path.isfile(spec.origin):
            candidates.append(spec.origin)
    except Exception:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(here, "libwlc.so"))
    candidates.extend(sorted(glob.glob(os.path.join(here, "_libwlc*.so"))))
    # Compiled in-place by setuptools (cpython tag / platform tag).
    candidates.extend(sorted(glob.glob(os.path.join(here, "_libwlc*.abi3.so"))))
    seen = set()
    for path in candidates:
        if not path or path in seen or not os.path.isfile(path):
            continue
        seen.add(path)
        lib = ctypes.CDLL(path)
        lib.wlc_encode.argtypes = [P, I, I, I, D, D, D, D, I, P, Lg, P, D]
        lib.wlc_encode.restype = Lg
        lib.wlc_decode.argtypes = [P, Lg, P]
        lib.wlc_shape.argtypes = [P, P]
        lib.wlc_encode_f.argtypes = [P, I, I, I, D, D, D, D, I, P, Lg, P, D]
        lib.wlc_encode_f.restype = Lg
        lib.wlc_decode_f.argtypes = [P, Lg, P]
        lib.wlc_info.argtypes = [P, P, P]
        lib.wlc_alloc.restype = ctypes.c_void_p
        lib.wlc_alloc.argtypes = [Lg]
        lib.wlc_free.argtypes = [P]
        lib.wlc_setfk.argtypes = [I]
        lib.wlc_setfk(FK)
        return lib
    raise ImportError(
        "tnvol native library not found. Install the package with pip "
        "(a C compiler and AVX2-capable x86-64 Linux are required) "
        "or build libwlc.so next to codec.py."
    )


lib = _load_lib()


def ptr(a):
    if isinstance(a, (bytearray, memoryview)):
        return ctypes.cast(
            ctypes.addressof((ctypes.c_char * len(a)).from_buffer(a)),
            ctypes.c_void_p,
        )
    if isinstance(a, int):
        return ctypes.c_void_p(a)
    return a.ctypes.data_as(ctypes.c_void_p)


def _kt(linf):
    return KT if not linf else KTL


# --- numpy-free .npy I/O for the CLI ---
def _hdr(t):
    d = t[t.index("'descr'"):]
    descr = d.split("'")[3]
    fo = "True" in t.split("'fortran_order'")[1].split(",")[0]
    shp = t.split("'shape'")[1]
    shp = shp[shp.index("(") + 1:shp.index(")")]
    return {
        "descr": descr,
        "fortran_order": fo,
        "shape": tuple(int(x) for x in shp.replace(" ", "").split(",") if x),
    }


def npy_read(path):
    import mmap
    with open(path, "rb") as f:
        m = f.read(8)
        if m[:6] != b"\x93NUMPY":
            raise ValueError("not a .npy file: %s" % path)
        hl = struct.unpack("<H", f.read(2))[0] if m[6] == 1 else struct.unpack("<I", f.read(4))[0]
        h = _hdr(f.read(hl).decode("latin1"))
        if h["fortran_order"]:
            raise ValueError("Fortran-order .npy is not supported; save as C-order")
        sh = tuple(h["shape"])
        n = 1
        for s_ in sh:
            n *= s_
        try:
            it = {"<u2": 2, "<f4": 4}[h["descr"]]
        except KeyError:
            raise ValueError("expects little-endian float32 or uint16, got %s" % h["descr"])
        off = f.tell()
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_COPY)
    return h["descr"], sh, memoryview(mm)[off:off + n * it]


def npy_write(path, descr, sh, buf):
    hd = "{'descr': '%s', 'fortran_order': False, 'shape': (%s), }" % (
        descr,
        ", ".join(map(str, sh)) + ("," if len(sh) == 1 else ""),
    )
    hd += " " * (64 - (10 + len(hd) + 1) % 64) + "\n"
    with open(path, "wb") as f:
        f.write(b"\x93NUMPY\x01\x00" + struct.pack("<H", len(hd)) + hd.encode("latin1"))
        f.write(buf)


def encode_raw(descr, sh, buf, psnr_t=None, linf=None):
    if len(sh) != 3:
        raise ValueError("volume must be 3-D, got shape %s" % (sh,))
    nz, ny, nx = sh
    n = nz * ny * nx
    st = (ctypes.c_double * 8)()
    if descr == "<f4":
        cap = n * 5 + (1 << 20)
        fn = lib.wlc_encode_f
    elif descr == "<u2":
        cap = n * 3 + (1 << 20)
        fn = lib.wlc_encode
    else:
        raise ValueError("expects uint16 or float32, got %s" % descr)
    out = lib.wlc_alloc(cap)
    if not out:
        raise MemoryError("wlc_alloc failed")
    try:
        k = fn(
            ptr(buf), nz, ny, nx,
            float(psnr_t or 0), float(linf or 0),
            _kt(linf), RND, MAXL,
            ctypes.c_void_p(out), cap, st, LAM,
        )
        if k <= 0:
            raise RuntimeError("encode failed (status %s)" % k)
        return bytes((ctypes.c_char * k).from_address(out))
    finally:
        lib.wlc_free(ctypes.c_void_p(out))


def decode_raw(blob):
    b = bytearray(blob)
    sh = (ctypes.c_int * 5)()
    tf = (ctypes.c_float * 1)()
    lib.wlc_info(ptr(b), sh, tf)
    n = sh[0] * sh[1] * sh[2]
    if n <= 0:
        raise ValueError("invalid stream header")
    it = 4 if sh[4] == 1 else 2
    r = lib.wlc_alloc(n * it)
    if not r:
        raise MemoryError("wlc_alloc failed")
    try:
        (lib.wlc_decode_f if sh[4] == 1 else lib.wlc_decode)(ptr(b), len(b), ctypes.c_void_p(r))
        return (
            "<f4" if it == 4 else "<u2",
            (sh[0], sh[1], sh[2]),
            bytes((ctypes.c_char * (n * it)).from_address(r)),
        )
    finally:
        lib.wlc_free(ctypes.c_void_p(r))


def encode(v, psnr_t=None, linf=None, stats=None):
    import numpy as np
    if np.asarray(v).dtype.kind == "f":
        v = np.ascontiguousarray(v, dtype=np.float32)
        nz, ny, nx = v.shape
        out = np.zeros(v.size * 5 + (1 << 20), np.uint8)
        st = np.zeros(8)
        n = lib.wlc_encode_f(
            ptr(v), nz, ny, nx,
            float(psnr_t or 0), float(linf or 0),
            _kt(linf), RND, MAXL,
            ptr(out), out.size, ptr(st), LAM,
        )
        if n <= 0:
            raise RuntimeError("encode failed (status %s)" % n)
        if stats is not None:
            stats[:] = st
        return out[:n].tobytes()
    v = np.ascontiguousarray(v, dtype=np.uint16)
    nz, ny, nx = v.shape
    out = np.zeros(v.size * 3 + (1 << 20), np.uint8)
    st = np.zeros(8)
    n = lib.wlc_encode(
        ptr(v), nz, ny, nx,
        float(psnr_t or 0), float(linf or 0),
        _kt(linf), RND, MAXL,
        ptr(out), out.size, ptr(st), LAM,
    )
    if n <= 0:
        raise RuntimeError("encode failed (status %s)" % n)
    if stats is not None:
        stats[:] = st
    return out[:n].tobytes()


def decode(blob):
    import numpy as np
    b = np.frombuffer(blob, np.uint8)
    sh = np.zeros(5, np.int32)
    tf = np.zeros(1, np.float32)
    lib.wlc_info(ptr(b), ptr(sh), ptr(tf))
    if sh[0] <= 0 or sh[1] <= 0 or sh[2] <= 0:
        raise ValueError("invalid stream header")
    if sh[4] == 1:
        r = np.empty((sh[0], sh[1], sh[2]), np.float32)
        lib.wlc_decode_f(ptr(b), b.size, ptr(r))
        return r, float(tf[0])
    r = np.empty((sh[0], sh[1], sh[2]), np.uint16)
    lib.wlc_decode(ptr(b), b.size, ptr(r))
    return r, int(sh[3])


def _par(jobs, fn):
    pids = []
    files = []
    try:
        for job in jobs:
            fd, path = tempfile.mkstemp(prefix="wlcmt")
            os.close(fd)
            files.append(path)
            pid = os.fork()
            if pid == 0:
                try:
                    open(path, "wb").write(fn(job))
                finally:
                    os._exit(0)
            pids.append(pid)
        for p in pids:
            os.waitpid(p, 0)
        out = []
        for path in files:
            with open(path, "rb") as f:
                out.append(f.read())
        return out
    finally:
        for path in files:
            try:
                os.remove(path)
            except OSError:
                pass


def encode_mt(descr, sh, buf, K, psnr_t=None, linf=None):
    import math
    import numpy as np
    nz, ny, nx = sh
    it = 4 if descr == "<f4" else 2
    dt = np.float32 if it == 4 else np.uint16
    v = np.frombuffer(buf, dt).reshape(sh)
    R = float(v.max()) - float(v.min())
    R = R if R > 0 else 1.0
    K = max(1, min(K, nz // 8))
    cuts = [round(i * nz / K) for i in range(K + 1)]
    jobs = []
    for i in range(K):
        sl = v[cuts[i]:cuts[i + 1]]
        Rs = float(sl.max()) - float(sl.min())
        Rs = Rs if Rs > 0 else 1.0
        p = max(1.0, psnr_t + 20 * math.log10(Rs / R)) if psnr_t else None
        jobs.append((cuts[i], cuts[i + 1], p))

    def fn(job):
        z0, z1, p = job
        m = memoryview(buf)[z0 * ny * nx * it:z1 * ny * nx * it]
        return bytes(encode_raw(descr, (z1 - z0, ny, nx), m, psnr_t=p, linf=linf))

    blobs = _par(jobs, fn)
    return b"WLCM" + struct.pack("<B", K) + b"".join(struct.pack("<Q", len(b)) for b in blobs) + b"".join(blobs)


def decode_mt(blob):
    K = blob[4]
    L = struct.unpack("<%dQ" % K, blob[5:5 + 8 * K])
    o = 5 + 8 * K
    parts = []
    for l in L:
        parts.append(blob[o:o + l])
        o += l
    res = _par(parts, lambda b: decode_raw(b)[2])
    shs = []
    for b in parts:
        h = (ctypes.c_int * 5)()
        tf = (ctypes.c_float * 1)()
        bb = bytearray(b[:4096] if len(b) > 4096 else b)
        lib.wlc_info(ptr(bb), h, tf)
        shs.append((h[0], h[1], h[2], h[4]))
    nz = sum(s[0] for s in shs)
    descr = "<f4" if shs[0][3] == 1 else "<u2"
    return descr, (nz, shs[0][1], shs[0][2]), b"".join(res)


def _require_target(psnr, linf):
    if (psnr is None or psnr == 0) and (linf is None or linf == 0):
        raise ValueError("specify psnr= (dB) or linf= (absolute L-infinity bound)")


def compress(volume, psnr=None, linf=None, threads=1):
    """Compress a 3-D numpy array (float32 or uint16).

    Provide either ``psnr`` (target peak-signal-to-noise ratio in dB) or
    ``linf`` (pointwise absolute error bound). Returns a compact bytes object.
    ``threads`` > 1 splits the volume into independent z-slabs (Linux only).
    """
    import numpy as np
    _require_target(psnr, linf)
    arr = np.asarray(volume)
    if arr.ndim != 3:
        raise ValueError("volume must be 3-D, got shape %s" % (arr.shape,))
    threads = int(threads or 1)
    if threads > 1:
        if arr.dtype.kind == "f":
            arr = np.ascontiguousarray(arr, dtype=np.float32)
            descr = "<f4"
        elif arr.dtype == np.uint16 or arr.dtype.kind in "iu":
            arr = np.ascontiguousarray(arr, dtype=np.uint16)
            descr = "<u2"
        else:
            raise TypeError("expects float32 or uint16, got %s" % arr.dtype)
        return encode_mt(descr, arr.shape, memoryview(arr), threads, psnr_t=psnr, linf=linf)
    return encode(arr, psnr_t=psnr, linf=linf)


def decompress(blob):
    """Decode a tnvol stream. Returns ``(array, tau)`` where ``tau`` is the
    L-infinity bound stored in the header (guaranteed ``|x - xhat| <= tau``
    at every voxel for a well-formed stream).
    """
    if not blob:
        raise ValueError("empty stream")
    if blob[:4] == b"WLCM":
        import numpy as np
        descr, sh, raw = decode_mt(blob)
        dt = np.float32 if descr == "<f4" else np.uint16
        rec = np.frombuffer(raw, dt).reshape(sh).copy()
        # Per-slab bounds can differ; report the max by re-reading headers.
        K = blob[4]
        L = struct.unpack("<%dQ" % K, blob[5:5 + 8 * K])
        o = 5 + 8 * K
        tau = 0.0
        for n in L:
            part = blob[o:o + n]
            o += n
            tau = max(tau, info(part)["tau"])
        return rec, tau
    return decode(blob)


def info(blob):
    """Return ``{shape, dtype, tau, multi}`` for a stream, without a full decode."""
    if blob[:4] == b"WLCM":
        K = blob[4]
        L = struct.unpack("<%dQ" % K, blob[5:5 + 8 * K])
        o = 5 + 8 * K
        first = None
        tau = 0.0
        nz = 0
        for n in L:
            part = blob[o:o + n]
            o += n
            meta = info(part)
            if first is None:
                first = meta
            nz += meta["shape"][0]
            tau = max(tau, meta["tau"])
        return {"shape": (nz, first["shape"][1], first["shape"][2]), "dtype": first["dtype"], "tau": tau, "multi": True}
    b = bytearray(blob[:4096] if len(blob) > 4096 else blob)
    sh = (ctypes.c_int * 5)()
    tf = (ctypes.c_float * 1)()
    lib.wlc_info(ptr(b), sh, tf)
    is_float = sh[4] == 1
    tau = float(tf[0]) if is_float else int(sh[3])
    return {
        "shape": (int(sh[0]), int(sh[1]), int(sh[2])),
        "dtype": "float32" if is_float else "uint16",
        "tau": tau,
        "multi": False,
    }


# Backwards-compatible aliases used in earlier snapshots.
def encode_volume(v, psnr_t=None, linf=None, stats=None):
    return encode(v, psnr_t=psnr_t, linf=linf, stats=stats)


def decode_volume(blob):
    return decode(blob)
