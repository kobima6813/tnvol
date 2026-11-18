# tnvol

**tnvol** is an error-bounded compressor for 3-D scientific volumes (seismic,
climate, turbulence, medical CT, and similar dense grids). It writes a compact
stream from a C-order `float32` or `uint16` array and reconstructs a volume
whose pointwise error stays within a bound stored in the header.

This package is the open-source core (frozen r10 codec). It is not a drop-in
replacement for every scientific compressor: please read the benchmark section,
including the cases where tnvol loses.

## What is in this package

- 3-D CDF 9/7 wavelet transform with packet (orientation) selection
- Dead-zone quantisation with rate-distortion optimisation
- Context-mixing binary range coder
- Wiener post-filter and a sparse L-infinity correction layer
- Palette / predictive mode for volumes with few distinct values

**Not in this release** (reserved for a commercial edition): brick-wise / ROI
random access, and a prediction mode specialised for sparse cloud fields.

## Requirements

- Linux on x86-64 with AVX2
- A C compiler (`gcc` or `clang`) for a source install
- Python 3.9+ and NumPy

## Install

```bash
pip install .
```

This compiles the native codec. The build uses `-DNPAIR=16 -mavx2 -mfma
-ffp-contract=off`, which is the configuration that matches the numbers below.
Do not change those flags if you want to reproduce the benches.

Nothing in this tree is published automatically. A public repository, if any,
will be announced separately.

## Quick start

### Command line

```bash
# Target a PSNR (dB)
tnvol compress volume.npy volume.tnv --psnr 45

# Or a pointwise absolute bound: |x - xhat| <= 0.01 at every voxel
tnvol compress volume.npy volume.tnv --linf 0.01

tnvol decompress volume.tnv reconstructed.npy
tnvol info volume.tnv
```

`volume.npy` must be a C-order 3-D array of little-endian `float32` or `uint16`.

### Python

```python
import numpy as np
import tnvol

vol = np.load("volume.npy")          # shape (Z, Y, X), float32 or uint16
blob = tnvol.compress(vol, psnr=45)  # or linf=0.01
rec, tau = tnvol.decompress(blob)

assert rec.shape == vol.shape
assert np.max(np.abs(vol.astype(np.float64) - rec.astype(np.float64))) <= tau
```

In `--linf` / `linf=` mode the encoder guarantees the bound you asked for. In
`--psnr` mode the header stores the achieved maximum absolute error (the
correction layer uses that value); the reconstructed volume stays within it.

Optional `--threads N` (or `threads=N`) splits the volume into independent
z-slabs. That path is Linux-only and changes the stream layout (`WLCM`
container). The default single-thread stream is the one used in the benches.

## Benchmarks

**Negative BD-rate means tnvol is smaller (better).** Figures are from the
frozen r10 codec (18 November 2026). Held-out, v3 and v4 tables were first
measured on the r9 wavelet/palette-less path and then checked: r10 produced
byte-identical streams on all 60 of those volumes, so the BD-rates carry over.
Seismic and gate numbers are from r10 itself.

BD-rate is computed on *measured* PSNR, using the Pareto front of each codec's
runs, not on the requested target. That matters for TTHRESH, which often
overshoots the requested PSNR (on F3 a request of 45 dB lands near 53 dB).

### Mean BD-rate by set

| Set | vs SPERR @PSNR | vs SZ3 @PSNR | vs TTHRESH @PSNR | vs SPERR @L-inf | vs SZ3 @L-inf |
|---|---:|---:|---:|---:|---:|
| Held-out (21 vols) | −23.1% | −38.3% | −42.2% | −34.5% | −25.5% |
| v3 (16 vols) | −20.7% | −48.6% | −44.6% | −36.7% | −37.6% |
| v4 (8 vols) | −30.6% | −53.6% | — | −42.2% | −40.3% |

Gate set (15 volumes: 6 Medical Segmentation Decathlon CTs + 9 SDRBench
fields): every volume passed the internal size/time gate. BD-rate vs SPERR at
PSNR is about −21% to −26% on the medical volumes and −20% to −41% on the
SDRBench volumes; vs TTHRESH at PSNR it is about −19% to −53% (medical) and
−27% to −79% (SDRBench). See `benches/RESULTS.md` for the per-volume table.

### Seismic volumes (BD-rate @PSNR)

| Volume | vs SPERR | vs SZ3 | vs TTHRESH |
|---|---:|---:|---:|
| F3 | −10.3% | −20.2% | **+14.4%** |
| Overthrust | −61.3% | −37.6% | −81.8% |
| Parihaka | −33.6% | −30.6% | −39.7% |
| Penobscot | −33.2% | −34.6% | −36.5% |
| Salt | −74.8% | −10.6% | −92.0% |

On F3, TTHRESH wins on the PSNR curve as a whole. A single low-rate point is
kinder to tnvol but is not the headline: at a matched PSNR of 45 dB, r10 is
0.7878 bits/voxel, which is 17% smaller than TTHRESH (interpolated to the same
quality; TTHRESH's raw “PSNR 45” run is not comparable because it achieved
52.9 dB) and 43% smaller than SPERR at the same requested point (1.3773
bits/voxel). The lead over both shrinks at high quality, which is why the
F3 BD-rate vs TTHRESH is still +14.4%.

L-inf BD-rates and CPU times are in `benches/RESULTS.md`. L-inf bounds were
satisfied on all reported r10 points. Encode+decode time is typically around
0.6–1.5× SPERR, slower than SZ3, and much faster than TTHRESH.

### Where tnvol loses

These are real losses, not footnotes:

- **F3 vs TTHRESH:** +14.4% BD-rate at PSNR (TTHRESH is smaller on the curve).
- **Sparse cloud / hydrometeor fields vs SZ3.** SZ3's Lorenzo predictor is a
  better fit for these near-sparse fields. Examples of BD-rate vs SZ3:
  - `h_hurr_CLOUDlog`: +67.1% @PSNR, +70.3% @L-inf
  - `v_hurr_QICElog`: +63.6% @PSNR, +73.8% @L-inf
  - `w_hurr_QCLOUDlog`: +76.3% @PSNR, +92.2% @L-inf
  - also `v_hurr_QRAINlog`, `v_hurr_QGRAUPlog`, `h_hurr_QSNOWlog` at L-inf
- **QMCPack (`h_qmc`) vs TTHRESH:** +28.0% BD-rate at PSNR.

On those same cloud fields tnvol is still smaller than SPERR. The open-source
core does not include the separate cloud prediction mode mentioned above.

Per-volume tables for every set, including ZFP and MGARD, are in
`benches/RESULTS.md`.

## How to reproduce

This package does not ship the scientific volumes or the competitor binaries
(we do not have redistribution rights for all of the datasets). To rerun the
comparison:

1. Build this package with the default compile flags (see Install). Do not
   change `NPAIR`, rounding, or the default `WLC_*` encoder knobs.
2. Obtain the volumes:
   - [SDRBench](https://sdrbench.github.io/) fields (Miranda, NYX, Hurricane
     Isabel, SCALE-LETKF, CESM, QMCPack, and the held-out / v3 / v4 splits of
     the same sources)
   - [Medical Segmentation Decathlon](http://medicaldecathlon.com/) CT volumes
     (liver, spleen, prostate), stored as 3-D `float32` or `uint16` `.npy`
   - Public seismic surveys: F3 Netherlands, Parihaka, Penobscot, SEG/EAGE
     Overthrust, SEAM salt, as 3-D `float32` grids in C-order
3. Compress with `tnvol compress --psnr {45,55,65,75,85}` and with `--linf`
   at the same absolute taus used for the bounded competitors (taus scale with
   each volume's range).
4. Run SZ3, SPERR, TTHRESH, and optionally ZFP and MGARD, at the same requested
   targets. Record bits/voxel, *measured* PSNR, measured L-inf, and CPU time.
5. Compute Bjøntegaard BD-rate on measured PSNR (Pareto front of each codec).
   For a fair point comparison, interpolate a competitor that overshoots the
   requested PSNR onto the same measured quality; do not compare raw
   “requested 45 dB” bitrates when the achieved PSNR differs by several dB.

A compact record of the numbers we measured is in `benches/RESULTS.md`.

## Limitations

- Supported platform today: Linux x86-64 with AVX2. There is no scalar fallback.
- Inputs must be finite: NaN and Inf are not supported.
- Only 3-D C-order `float32` and `uint16` are accepted.
- Multi-slab `--threads` uses `fork()` and is not a bit-identical substitute
  for the single-thread stream used in the benches.

## License

Copyright (c) 2026 Yackov Mozgori.

tnvol is free software under the **GNU Affero General Public License v3.0**
only (`LICENSE`). If you cannot comply with the AGPL — for example if you want
to use tnvol in a proprietary product or a hosted service without releasing
your own source — you can buy a commercial license. See `COMMERCIAL.md`.
