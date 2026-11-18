# SPDX-License-Identifier: AGPL-3.0-only
"""Frozen r10 bit-identity: encodings must match the published r10 streams.

Hashes were taken from the frozen r10 libwlc.so (18 November 2026) on the
volumes used by test_roundtrip.py plus a few random float32/uint16 cases.
On x86-64, if both the portable and AVX2 libraries are present, they must
also match each other byte-for-byte (encode stream and decoded array).
"""
from __future__ import annotations

import hashlib
import os
import sys
import unittest

import numpy as np

import tnvol
from tnvol import codec

# sha256 of r10 streams produced with the default encoder knobs
# (RND=0.38, MAXL=6, LAM=0.08, FK=6, kt=KTL for linf / KT for psnr).
GOLDEN = {
    "smooth_linf": "d603188de2891d84de411c6312f374c8b7bde7955bc784b968b903dd9e810aea",
    "smooth_psnr45": "b6c1ed2b25c5157d0a8704301be10c1b479215cf335b1a8df43fe0d0caf819f0",
    "u16_linf": "a0405f63c6d1b0c92980afe0efa4d3c8509d6d6990d7e66b6f760bee1be55dd6",
    "blocky_linf": "227d5df12c6241406855231a8cf3f364cbcd1aa590785785fc67fc9b5ac8d2f1",
    "smooth16_linf": "13b43de28deeb9e754c433c99c40d68b43b697260883dab6e5a8c3df824b2199",
    "smooth24_linf": "803066272cb995586c3288bdcad40e41dcebf354baa552f3c1b88e1b3473a61a",
    "randf_20x20x20_linf": "fbcaffddaab0371ade1ad17ed87a650f3e17a3640d2f0e504634bd39139cb9fa",
    "randf_20x20x20_psnr": "669b23dfef9c4f26e5df6148c18487a66486ef75f1d64db28ad0ec498486da1a",
    "randu_20x20x20_linf": "d30bdf178a28987ae82791406e6b00d943c93dfc52c9a01c04d814dd55113ce1",
    "randf_40x16x24_linf": "330f2991fc040a1674c5fb22681adb69d03ab05b80ff9bd759c9ba8c02d81476",
    "randf_40x16x24_psnr": "41be7b2d7709081bc040f73e55bfe24e0dab03cd819603406b7b8d19e3c8792a",
    "randu_40x16x24_linf": "a6fa06bcacc93433f98f926828c2d371631cd14c6ed82897d1ac39b4dd591505",
    "randf_12x48x36_linf": "050ccf33d4943c52600d19019beab4cab105918dbde17bdf7f4ad7ccd35110ee",
    "randf_12x48x36_psnr": "03331c0caa24e14bfebb96ea287de17f6478b9148896022ea6ad868659b3b529",
    "randu_12x48x36_linf": "2de40cbbf2c901525ccf2150975abdf54f04f86cc0344778d571945d6f2814f1",
    "randf_8x8x8_linf": "77febaa19632d72f5b7fcf0f96f8c280afe1e73f5bca8eca7e69d7a8d39717f3",
    "randf_8x8x8_psnr": "461eabffc924b69db50a99d8a7bebd9a2c70e5cbf80f8c4dd5968e59a8fbce6f",
    "randu_8x8x8_linf": "4d3d8267a0cae0f489bf7da7f112de9afa2a7c8a75b1be78bcb63308eeb565cb",
    "randf_33x31x17_linf": "304f4f72c766dcea319ed7eb49a3a6a3f21cc2d52c28af0036263a4cc90408ac",
    "randf_33x31x17_psnr": "1fd451d3c9849b835b238acb0ca92dd06e79e341243827293adb1ff32a8f7509",
    "randu_33x31x17_linf": "5b256335d1348ff67c6ef76611689b638fbf6933ee57d83c6bdc666254f91b92",
}


def _smooth(shape=(32, 40, 48), seed=0):
    rng = np.random.default_rng(seed)
    z, y, x = np.mgrid[0:shape[0], 0:shape[1], 0:shape[2]]
    vol = (np.sin(z / 5.0) + np.cos(y / 7.0) + np.sin(x / 6.0)).astype(np.float32)
    vol += 0.05 * rng.standard_normal(shape, dtype=np.float32)
    return vol


def _blocky(shape=(24, 32, 32)):
    vol = np.zeros(shape, dtype=np.float32)
    vol[5:15, 8:20, 8:20] = 3.5
    vol[10:18, 12:28, 4:16] = -1.25
    vol[0, 0, 0] = 0.5
    return vol


def _cases():
    cases = [
        ("smooth_linf", _smooth(), dict(linf=0.05)),
        ("smooth_psnr45", _smooth(), dict(psnr=45)),
        ("u16_linf", (1000 + 200 * _smooth((32, 32, 32))).clip(0, 65535).astype(np.uint16), dict(linf=4)),
        ("blocky_linf", _blocky(), dict(linf=0.01)),
        ("smooth16_linf", _smooth((16, 24, 24)), dict(linf=0.1)),
        ("smooth24_linf", _smooth((24, 24, 24)), dict(linf=0.02)),
    ]
    rng = np.random.default_rng(123)
    for shape in [(20, 20, 20), (40, 16, 24), (12, 48, 36), (8, 8, 8), (33, 31, 17)]:
        vol = rng.standard_normal(shape, dtype=np.float32)
        cases.append((f"randf_{shape[0]}x{shape[1]}x{shape[2]}_linf", vol, dict(linf=0.08)))
        cases.append((f"randf_{shape[0]}x{shape[1]}x{shape[2]}_psnr", vol, dict(psnr=50)))
        u16 = rng.integers(0, 65535, size=shape, dtype=np.uint16)
        cases.append((f"randu_{shape[0]}x{shape[1]}x{shape[2]}_linf", u16, dict(linf=8)))
    return cases


def _sibling_libs():
    here = os.path.dirname(os.path.abspath(codec.__file__))
    exts = {".so", ".dylib", ".dll"}
    port, avx = [], []
    for name in os.listdir(here):
        root, ext = os.path.splitext(name)
        if ext.lower() not in exts:
            continue
        path = os.path.join(here, name)
        if "avx2" in name.lower():
            avx.append(path)
        elif name.startswith("libwlc") or name.lower() in ("wlc.dll", "libwlc.dll"):
            port.append(path)
    return port, avx


class TestBitIdentity(unittest.TestCase):
    def test_frozen_r10_hashes(self):
        for name, vol, kw in _cases():
            blob = tnvol.compress(vol, **kw)
            got = hashlib.sha256(blob).hexdigest()
            self.assertEqual(got, GOLDEN[name], msg="r10 mismatch for %s" % name)
            rec, _ = tnvol.decompress(blob)
            self.assertEqual(rec.shape, vol.shape)
            self.assertEqual(rec.dtype, vol.dtype)

    def test_portable_matches_avx2_when_both_present(self):
        import ctypes

        port, avx = _sibling_libs()
        if not port or not avx:
            self.skipTest("both portable and AVX2 libraries are not installed")
        def load(path):
            lib = ctypes.CDLL(path)
            P, I, Lg, D = ctypes.c_void_p, ctypes.c_int, ctypes.c_long, ctypes.c_double
            lib.wlc_encode.argtypes = [P, I, I, I, D, D, D, D, I, P, Lg, P, D]
            lib.wlc_encode.restype = Lg
            lib.wlc_decode.argtypes = [P, Lg, P]
            lib.wlc_encode_f.argtypes = [P, I, I, I, D, D, D, D, I, P, Lg, P, D]
            lib.wlc_encode_f.restype = Lg
            lib.wlc_decode_f.argtypes = [P, Lg, P]
            lib.wlc_setfk.argtypes = [I]
            lib.wlc_setfk(6)
            return lib

        libs = {"port": load(port[0]), "avx2": load(avx[0])}

        def run(lib, vol, **kw):
            nz, ny, nx = vol.shape
            cap = max(vol.nbytes * 5 + (1 << 20), 2_000_000)
            out = np.zeros(cap, np.uint8)
            st = np.zeros(8)
            kt = float(os.environ.get("WLC_KTL", -1)) if kw.get("linf") else float(os.environ.get("WLC_KT", 0.0))
            rnd = float(os.environ.get("WLC_RND", 0.38))
            maxl = int(os.environ.get("WLC_MAXL", 6))
            lam = float(os.environ.get("WLC_LAM", 0.08))
            fn = lib.wlc_encode_f if vol.dtype == np.float32 else lib.wlc_encode
            n = fn(
                vol.ctypes.data_as(ctypes.c_void_p), nz, ny, nx,
                float(kw.get("psnr") or 0), float(kw.get("linf") or 0),
                kt, rnd, maxl, out.ctypes.data_as(ctypes.c_void_p), cap,
                st.ctypes.data_as(ctypes.c_void_p), lam,
            )
            self.assertGreater(n, 0)
            blob = bytes(out[:n])
            rec = np.empty(vol.shape, vol.dtype)
            dfn = lib.wlc_decode_f if vol.dtype == np.float32 else lib.wlc_decode
            buf = (ctypes.c_char * len(blob)).from_buffer_copy(blob)
            dfn(ctypes.cast(buf, ctypes.c_void_p), len(blob), rec.ctypes.data_as(ctypes.c_void_p))
            return blob, rec

        for name, vol, kw in _cases():
            b0, r0 = run(libs["port"], vol, **kw)
            b1, r1 = run(libs["avx2"], vol, **kw)
            self.assertEqual(b0, b1, msg="encode bytes differ for %s" % name)
            self.assertTrue(np.array_equal(r0, r1), msg="decode arrays differ for %s" % name)


if __name__ == "__main__":
    unittest.main()
