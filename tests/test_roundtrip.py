# SPDX-License-Identifier: AGPL-3.0-only
"""Round-trip and error-bound checks. Run after ``pip install .``."""
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import tnvol


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


class TestRoundtrip(unittest.TestCase):
    def test_linf_float32_bound(self):
        vol = _smooth()
        tau = 0.05
        blob = tnvol.compress(vol, linf=tau)
        rec, stored = tnvol.decompress(blob)
        self.assertEqual(rec.shape, vol.shape)
        self.assertEqual(rec.dtype, np.float32)
        err = float(np.max(np.abs(vol - rec)))
        self.assertLessEqual(err, tau + 1e-6)
        self.assertLessEqual(err, stored + 1e-6)
        self.assertGreater(len(blob), 16)

    def test_psnr_float32_bound(self):
        vol = _smooth()
        blob = tnvol.compress(vol, psnr=45)
        rec, stored = tnvol.decompress(blob)
        err = float(np.max(np.abs(vol - rec)))
        self.assertLessEqual(err, stored + 1e-5)
        # A 45 dB target should be in the right ballpark on this smooth field.
        rms = float(np.sqrt(np.mean((vol - rec) ** 2)))
        peak = float(vol.max() - vol.min())
        psnr = 20.0 * np.log10(peak / max(rms, 1e-30))
        self.assertGreater(psnr, 40.0)

    def test_linf_uint16_bound(self):
        vol = (1000 + 200 * _smooth((32, 32, 32))).clip(0, 65535).astype(np.uint16)
        tau = 4
        blob = tnvol.compress(vol, linf=tau)
        rec, stored = tnvol.decompress(blob)
        self.assertEqual(rec.dtype, np.uint16)
        err = int(np.max(np.abs(vol.astype(np.int32) - rec.astype(np.int32))))
        self.assertLessEqual(err, int(tau))
        self.assertLessEqual(err, int(stored))

    def test_palette_like_volume(self):
        vol = _blocky()
        tau = 0.01
        blob = tnvol.compress(vol, linf=tau)
        rec, stored = tnvol.decompress(blob)
        err = float(np.max(np.abs(vol - rec)))
        self.assertLessEqual(err, tau + 1e-6)
        self.assertLessEqual(err, stored + 1e-6)
        # A few distinct values should compress to a very small stream.
        self.assertLess(len(blob), vol.nbytes // 8)

    def test_info(self):
        vol = _smooth((16, 24, 24))
        blob = tnvol.compress(vol, linf=0.1)
        meta = tnvol.info(blob)
        self.assertEqual(meta["shape"], vol.shape)
        self.assertEqual(meta["dtype"], "float32")
        self.assertGreater(meta["tau"], 0)

    def test_cli_roundtrip(self):
        vol = _smooth((24, 24, 24))
        with tempfile.TemporaryDirectory() as td:
            src = os.path.join(td, "in.npy")
            packed = os.path.join(td, "out.tnv")
            rec_path = os.path.join(td, "rec.npy")
            np.save(src, vol)
            subprocess.check_call(
                [sys.executable, "-m", "tnvol.cli", "compress", src, packed, "--linf", "0.02"]
            )
            subprocess.check_call(
                [sys.executable, "-m", "tnvol.cli", "decompress", packed, rec_path]
            )
            rec = np.load(rec_path)
            err = float(np.max(np.abs(vol - rec)))
            self.assertLessEqual(err, 0.02 + 1e-6)


if __name__ == "__main__":
    unittest.main()
