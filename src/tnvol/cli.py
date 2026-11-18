# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Yackov Mozgori
"""Command-line interface: ``tnvol compress`` / ``tnvol decompress``."""
from __future__ import annotations

import argparse
import os
import sys


def _build_parser():
    p = argparse.ArgumentParser(
        prog="tnvol",
        description="Error-bounded compressor for 3-D scientific volumes (float32 or uint16 .npy).",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("compress", help="compress a C-order .npy volume")
    c.add_argument("volume", nargs="?", help="input .npy (C-order float32 or uint16)")
    c.add_argument("out", nargs="?", help="output stream path")
    c.add_argument("--volume", dest="volume_opt", help=argparse.SUPPRESS)
    c.add_argument("--out", dest="out_opt", help=argparse.SUPPRESS)
    c.add_argument("--psnr", type=float, default=None, help="target PSNR in dB")
    c.add_argument("--linf", type=float, default=None, help="pointwise absolute error bound")
    c.add_argument("--threads", type=int, default=None, help="independent z-slabs (Linux); default 1")
    c.add_argument("--effort", type=int, default=None, help="encoder effort (2 = post-filter-aware requantisation)")

    d = sub.add_parser("decompress", help="decompress a stream to .npy")
    d.add_argument("packed", nargs="?", help="input stream")
    d.add_argument("out", nargs="?", help="output .npy path")
    d.add_argument("--packed", dest="packed_opt", help=argparse.SUPPRESS)
    d.add_argument("--out", dest="out_opt", help=argparse.SUPPRESS)

    i = sub.add_parser("info", help="print shape, dtype and stored error bound")
    i.add_argument("packed", nargs="?", help="input stream")
    i.add_argument("--packed", dest="packed_opt", help=argparse.SUPPRESS)
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    p = _build_parser()
    args = p.parse_args(argv)

    # Lazy import so ``tnvol --help`` works even if the native lib failed to build.
    if args.cmd == "compress":
        volume = args.volume or args.volume_opt
        out = args.out or args.out_opt
        if not volume or not out:
            p.error("compress requires an input .npy and an output path")
        if args.psnr is None and args.linf is None:
            p.error("specify --psnr or --linf")
        if args.effort is not None:
            os.environ["WLC_EFFORT"] = str(args.effort)
        from .codec import npy_read, encode_raw, encode_mt
        descr, sh, buf = npy_read(volume)
        threads = args.threads
        if threads is None:
            threads = int(os.environ.get("WLC_THREADS", 1))
        if threads > 1:
            blob = encode_mt(descr, sh, buf, threads, psnr_t=args.psnr, linf=args.linf)
        else:
            blob = encode_raw(descr, sh, buf, psnr_t=args.psnr, linf=args.linf)
        with open(out, "wb") as f:
            f.write(blob)
        return 0

    if args.cmd == "decompress":
        packed = args.packed or args.packed_opt
        out = args.out or args.out_opt
        if not packed or not out:
            p.error("decompress requires an input stream and an output .npy path")
        from .codec import npy_write, decode_raw, decode_mt
        with open(packed, "rb") as f:
            blob = f.read()
        if blob[:4] == b"WLCM":
            descr, sh, raw = decode_mt(blob)
        else:
            descr, sh, raw = decode_raw(blob)
        npy_write(out, descr, sh, raw)
        return 0

    if args.cmd == "info":
        packed = args.packed or getattr(args, "packed_opt", None)
        if not packed:
            p.error("info requires an input stream")
        from .codec import info
        with open(packed, "rb") as f:
            blob = f.read()
        meta = info(blob)
        print("shape=%s  dtype=%s  tau=%s  multi=%s" % (meta["shape"], meta["dtype"], meta["tau"], meta["multi"]))
        return 0

    p.error("unknown command")
    return 2


if __name__ == "__main__":
    sys.exit(main())
