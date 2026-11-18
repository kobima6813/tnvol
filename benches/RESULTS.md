# tnvol r10 benchmark numbers

Measured on the frozen r10 codec, 18 November 2026. Held-out, v3 and v4
BD-rates were first measured on r9; r10 produced byte-identical streams on all
60 of those volumes, so the rates below are the r10 figures as well. Seismic
and gate tables were measured on r10 directly.

**Negative BD-rate means tnvol is smaller (better).** BD-rate at PSNR uses
measured PSNR and the Pareto front of each codec's runs. BD-rate at L-inf uses
the same requested absolute tolerance for every bounded codec. "ours bound ok"
means the reconstructed volume satisfied `|err| <= tau` at every voxel.

CPU columns are mean encode+decode seconds on the machine used for the
original run; they are for orientation only.

## 1. Gate set (medical CT + SDRBench)

| volume | BD@PSNR vs SPERR | BD@L-inf vs SPERR | BD@L-inf vs SZ3 | BD@PSNR vs TTHRESH | CPU ours / SPERR |
|---|---:|---:|---:|---:|---:|
| la_007 | −25.4% | −28.7% | n/a | −28.5% | 0.65× |
| la_019 | −25.3% | −32.3% | n/a | −19.2% | 0.63× |
| la_023 | −23.6% | −27.5% | n/a | −26.1% | 0.67× |
| spleen_19 | −21.8% | −24.9% | n/a | −50.4% | 0.96× |
| spleen_31 | −23.4% | −27.4% | n/a | −48.8% | 0.97× |
| spleen_52 | −26.2% | −27.1% | n/a | −52.7% | 0.90× |
| miranda_density | −36.4% | −34.9% | −44.3% | −26.8% | 0.61× |
| miranda_pressure | −36.3% | −39.0% | −43.1% | −43.6% | 0.66× |
| miranda_velocityx | −36.6% | −41.5% | −46.6% | −35.9% | 0.62× |
| hurr_P | −40.9% | −35.7% | −36.3% | −74.9% | 0.92× |
| hurr_TC | −24.3% | −31.9% | −34.8% | −72.1% | 1.12× |
| hurr_U | −28.2% | −26.7% | −41.9% | −53.8% | 1.01× |
| nyx_temperature | −26.2% | −28.6% | −22.8% | −78.5% | 0.92× |
| nyx_velocityx | −28.6% | −32.2% | −18.7% | −67.4% | 0.92× |
| nyx_baryon_log10 | −20.5% | −20.1% | −15.5% | −41.9% | 1.04× |

All 15 volumes met the size and time criteria used in that run. SZ3 L-inf
BD-rate was not recorded for the six medical volumes in this table.

## 2. Held-out volumes

### BD-rate @PSNR

| volume | SPERR | SZ3 | ZFP | TTHRESH | MGARD |
|---|---:|---:|---:|---:|---:|
| h_la_011 | −19.9% | −37.2% | −31.2% | −24.2% | −56.4% |
| h_la_016 | −22.3% | −51.5% | −36.8% | −10.1% | −63.7% |
| h_spleen_6 | −13.2% | −36.0% | −29.0% | −28.2% | −46.3% |
| h_spleen_40 | −15.2% | −38.5% | −35.7% | −41.7% | −57.4% |
| h_prostate_00 | −12.6% | −17.6% | −30.8% | −26.5% | −35.1% |
| h_prostate_04 | −9.6% | −17.8% | −31.5% | −25.8% | −37.9% |
| h_hurr_W | −12.8% | −60.6% | −64.0% | −47.4% | −78.2% |
| h_hurr_CLOUDlog | −49.1% | **+67.1%** | −3.7% | −75.1% | +0.7% |
| h_hurr_QSNOWlog | −20.9% | −0.9% | −18.1% | −55.7% | −27.7% |
| h_hurr_PRECIP | −27.4% | −37.9% | −39.6% | −84.3% | −66.8% |
| h_nyx_vz | −15.1% | −36.1% | −37.9% | −57.0% | −57.3% |
| h_nyx_dm | −39.5% | −41.9% | −48.0% | −82.0% | −55.6% |
| h_mir_vz | −27.3% | −58.5% | −82.1% | −34.5% | −91.0% |
| h_mir_diff | −23.9% | −51.3% | −56.3% | −43.5% | −81.4% |
| h_scale_T | −19.0% | −89.5% | −89.9% | −63.0% | −95.0% |
| h_scale_U | −24.9% | −87.5% | −83.6% | −42.8% | −93.7% |
| h_scale_QRlog | −26.6% | −38.7% | −49.3% | −71.7% | −62.5% |
| h_qmc | −31.8% | −20.5% | −34.9% | **+28.0%** | −27.5% |
| h_basil_u | −25.5% | −78.6% | −94.2% | −16.4% | −96.9% |
| h_cesm_CLDHGH | −18.6% | −27.0% | −70.5% | — | — |
| h_cesm_FLDSC | −30.6% | −43.1% | −82.3% | — | — |
| **mean** | **−23.1%** | **−38.3%** | **−50.0%** | **−42.2%** | **−59.5%** |
| wins/n | 21/21 | 20/21 | 21/21 | 18/19 | 18/19 |

### BD-rate @L-inf (same requested tau)

| volume | SPERR | SZ3 | ZFP | MGARD | tnvol bound ok | CPU vs SPERR |
|---|---:|---:|---:|---:|---|---:|
| h_la_011 | −31.2% | −25.4% | −76.5% | −69.9% | yes | 0.70× |
| h_la_016 | −35.0% | −38.2% | −81.0% | −75.3% | yes | 0.76× |
| h_spleen_6 | −29.5% | −22.1% | −78.9% | −68.1% | yes | 1.05× |
| h_spleen_40 | −31.4% | −26.2% | −83.4% | −73.6% | yes | 0.98× |
| h_prostate_00 | −20.3% | −8.5% | −68.1% | −53.1% | yes | 1.29× |
| h_prostate_04 | −20.1% | −10.5% | −69.8% | −56.3% | yes | 1.37× |
| h_hurr_W | −27.0% | −45.0% | −88.0% | −82.4% | yes | 1.23× |
| h_hurr_CLOUDlog | −44.1% | **+70.3%** | −43.4% | −17.0% | yes | 2.14× |
| h_hurr_QSNOWlog | −22.9% | **+7.5%** | −55.1% | −41.2% | yes | 1.26× |
| h_hurr_PRECIP | −37.7% | −24.4% | −83.2% | −76.7% | yes | 1.18× |
| h_nyx_vz | −28.4% | −17.5% | −83.2% | −70.2% | yes | 1.18× |
| h_nyx_dm | −55.6% | −36.1% | −91.5% | −77.9% | yes | 1.18× |
| h_mir_vz | −43.1% | −40.8% | −93.6% | −91.4% | yes | 0.73× |
| h_mir_diff | −37.0% | −36.4% | −85.6% | −85.0% | yes | 0.81× |
| h_scale_T | −45.2% | −67.3% | −97.7% | −93.2% | yes | 0.95× |
| h_scale_U | −45.2% | −70.8% | −96.5% | −93.3% | yes | 1.07× |
| h_scale_QRlog | −31.8% | −21.3% | −80.1% | −68.7% | yes | 1.31× |
| h_qmc | −35.4% | −20.2% | −63.2% | −45.5% | yes | 0.84× |
| h_basil_u | −46.7% | −67.9% | −97.5% | −96.7% | yes | 0.86× |
| h_cesm_CLDHGH | −18.8% | −11.1% | −90.9% | — | yes | 1.12× |
| h_cesm_FLDSC | −38.8% | −22.7% | −96.7% | — | yes | 1.00× |
| **mean** | **−34.5%** | **−25.5%** | **−81.1%** | **−70.3%** |  |  |
| wins/n | 21/21 | 19/21 | 21/21 | 19/19 |  |  |

## 3. v3 volumes

### BD-rate @PSNR

| volume | SPERR | SZ3 | ZFP | TTHRESH | MGARD |
|---|---:|---:|---:|---:|---:|
| v_la_014 | −14.8% | −41.7% | −31.9% | −26.3% | −57.6% |
| v_la_021 | −14.2% | −42.8% | −33.4% | −21.1% | −61.1% |
| v_spleen_8 | −10.2% | −29.3% | −28.2% | −37.8% | −40.1% |
| v_spleen_44 | −8.3% | −41.6% | −31.6% | −58.4% | −51.8% |
| v_prostate_02 | −9.7% | −11.3% | −20.1% | −17.2% | −24.8% |
| v_scale_PRES | −30.5% | −91.8% | −98.2% | −73.2% | −98.0% |
| v_scale_V | −19.9% | −86.2% | −81.6% | −35.1% | −92.6% |
| v_scale_RH | −16.9% | −81.4% | −72.8% | −35.9% | −89.5% |
| v_scale_QVlog | −21.3% | −89.5% | −88.4% | −48.4% | −95.5% |
| v_hurr_QGRAUPlog | −21.5% | **+0.5%** | −16.3% | −54.1% | −26.5% |
| v_hurr_QICElog | −51.9% | **+63.6%** | −10.2% | −79.3% | +0.7% |
| v_hurr_QRAINlog | −19.2% | **+7.2%** | −12.4% | −47.7% | −22.4% |
| v_basil_v | −26.0% | −80.6% | −94.1% | −19.9% | −97.1% |
| v_cesm_T | −21.3% | −91.3% | −93.7% | −65.7% | −96.8% |
| v_cesm_U | −20.5% | −87.8% | −87.3% | −53.9% | −95.2% |
| v_cesm_RELHUM | −24.9% | −74.3% | −70.7% | −39.1% | −86.7% |
| **mean** | **−20.7%** | **−48.6%** | **−54.4%** | **−44.6%** | **−64.7%** |
| wins/n | 16/16 | 13/16 | 16/16 | 16/16 | 15/16 |

### BD-rate @L-inf

| volume | SPERR | SZ3 | ZFP | MGARD | tnvol bound ok | CPU vs SPERR |
|---|---:|---:|---:|---:|---|---:|
| v_la_014 | −27.5% | −29.1% | −76.8% | −70.5% | yes | 0.78× |
| v_la_021 | −28.6% | −30.5% | −78.3% | −72.9% | yes | 0.73× |
| v_spleen_8 | −25.0% | −20.7% | −75.3% | −63.9% | yes | 1.06× |
| v_spleen_44 | −29.2% | −33.6% | −80.7% | −71.6% | yes | 0.88× |
| v_prostate_02 | −16.0% | −4.7% | −58.9% | −43.6% | yes | 1.04× |
| v_scale_PRES | −66.6% | −79.9% | −99.5% | −97.0% | yes | 1.02× |
| v_scale_V | −42.6% | −69.8% | −95.8% | −92.6% | yes | 1.22× |
| v_scale_RH | −37.1% | −64.3% | −94.0% | −90.1% | yes | 1.30× |
| v_scale_QVlog | −49.6% | −71.4% | −97.7% | −94.5% | yes | 1.20× |
| v_hurr_QGRAUPlog | −22.9% | **+8.7%** | −53.2% | −39.7% | yes | 1.45× |
| v_hurr_QICElog | −46.8% | **+73.8%** | −49.4% | −20.9% | yes | 2.42× |
| v_hurr_QRAINlog | −21.7% | **+12.9%** | −50.1% | −37.2% | yes | 8.00× |
| v_basil_v | −47.2% | −70.7% | −97.4% | −96.9% | yes | 1.03× |
| v_cesm_T | −48.3% | −82.4% | −98.1% | −96.7% | yes | 1.27× |
| v_cesm_U | −42.9% | −77.4% | −96.8% | −95.5% | yes | 1.40× |
| v_cesm_RELHUM | −35.0% | −62.7% | −91.5% | −89.2% | yes | 1.50× |
| **mean** | **−36.7%** | **−37.6%** | **−80.9%** | **−73.3%** |  |  |
| wins/n | 16/16 | 13/16 | 16/16 | 16/16 |  |  |

The 8.00× CPU ratio on `v_hurr_QRAINlog` is an outlier (27 s vs 3.4 s SPERR on
that run). Other v3 volumes are near 1× SPERR.

## 4. v4 volumes

### BD-rate @PSNR

| volume | SPERR | SZ3 | ZFP | MGARD |
|---|---:|---:|---:|---:|
| w_scale_W | −16.0% | −76.9% | −75.7% | −88.0% |
| w_scale_QSlog | −41.1% | −51.6% | −62.4% | −73.7% |
| w_hurr_QCLOUDlog | −66.8% | **+76.3%** | −20.3% | −19.2% |
| w_cesm_Z3 | −33.5% | −93.5% | −97.7% | −97.1% |
| w_cesm_OMEGA | −18.8% | −61.6% | −57.4% | −77.1% |
| w_cesm_Q | −19.5% | −86.7% | −85.3% | −93.2% |
| w_cesm_CLOUD | −31.0% | −51.6% | −56.0% | −70.7% |
| w_cesm_V | −18.4% | −82.9% | −82.0% | −92.4% |
| **mean** | **−30.6%** | **−53.6%** | **−67.1%** | **−76.4%** |
| wins/n | 8/8 | 7/8 | 8/8 | 8/8 |

TTHRESH was not recorded for this set.

### BD-rate @L-inf

| volume | SPERR | SZ3 | ZFP | MGARD | tnvol bound ok | CPU vs SPERR |
|---|---:|---:|---:|---:|---|---:|
| w_scale_W | −38.3% | −58.8% | −94.6% | −89.0% | yes | 1.16× |
| w_scale_QSlog | −43.6% | −32.2% | −85.6% | −77.2% | yes | 1.37× |
| w_hurr_QCLOUDlog | −57.4% | **+92.2%** | −55.7% | −18.3% | yes | 2.50× |
| w_cesm_Z3 | −62.7% | −87.4% | −99.3% | −96.8% | yes | 1.08× |
| w_cesm_OMEGA | −26.9% | −48.1% | −87.7% | −82.1% | yes | 1.64× |
| w_cesm_Q | −38.8% | −76.2% | −96.0% | −93.5% | yes | 1.49× |
| w_cesm_CLOUD | −32.8% | −40.3% | −83.5% | −77.7% | yes | 1.57× |
| w_cesm_V | −36.7% | −72.1% | −94.8% | −93.4% | yes | 1.82× |
| **mean** | **−42.2%** | **−40.3%** | **−87.1%** | **−78.5%** |  |  |
| wins/n | 8/8 | 7/8 | 8/8 | 8/8 |  |  |

## 5. Seismic volumes (r10)

### Bits/voxel at requested targets

SPERR column is shown where it was recorded (PSNR mode). SZ3 is shown for
L-inf mode. "stream" is the r10 mode that won (wavelet vs palette).

| volume | mode | target | r10 bpv | stream | L-inf ok | SZ3 bpv | SPERR bpv |
|---|---|---:|---:|---|---|---:|---:|
| f3 | linf | 2e-05 | 10.9572 | wavelet | yes | 11.5827 | 11.4393 |
| f3 | linf | 0.0002 | 7.6088 | wavelet | yes | 8.2579 | 8.1197 |
| f3 | linf | 0.002 | 4.1646 | wavelet | yes | 4.9223 | 4.7953 |
| f3 | linf | 0.02 | 1.1890 | wavelet | yes | 1.6726 | 1.7070 |
| f3 | psnr | 45 | 0.7878 | wavelet |  |  | 1.3773 |
| f3 | psnr | 55 | 2.0698 | wavelet |  |  | 2.6160 |
| f3 | psnr | 65 | 3.8183 | wavelet |  |  | 4.1869 |
| f3 | psnr | 75 | 5.5328 | wavelet |  |  | 5.9735 |
| f3 | psnr | 85 | 7.2250 | wavelet |  |  | 7.7129 |
| overthrust | linf | 0.0373035 | 2.3974 | wavelet | yes | 2.4727 | 4.8759 |
| overthrust | linf | 0.373035 | 1.4121 | wavelet | yes | 1.5688 | 3.1345 |
| overthrust | linf | 3.73035 | 0.6456 | wavelet | yes | 0.8771 | 1.6349 |
| overthrust | linf | 37.3035 | 0.1831 | wavelet | yes | 0.3175 | 0.5897 |
| overthrust | psnr | 45 | 0.0432 | wavelet |  |  | 0.5038 |
| overthrust | psnr | 55 | 0.1569 | wavelet |  |  | 0.9120 |
| overthrust | psnr | 65 | 0.3732 | wavelet |  |  | 1.4556 |
| overthrust | psnr | 75 | 0.6883 | wavelet |  |  | 2.1257 |
| overthrust | psnr | 85 | 1.0747 | wavelet |  |  | 2.8904 |
| parihaka | linf | 0.00255 | 4.6970 | palette | yes | 6.8482 | 12.9431 |
| parihaka | linf | 0.0255 | 4.6970 | palette | yes | 6.9005 | 9.6215 |
| parihaka | linf | 0.255 | 4.6970 | palette | yes | 6.1782 | 6.2986 |
| parihaka | linf | 2.55 | 2.2075 | wavelet | yes | 2.8182 | 2.9726 |
| parihaka | psnr | 45 | 1.8134 | wavelet |  |  | 2.4862 |
| parihaka | psnr | 55 | 3.3871 | wavelet |  |  | 4.0264 |
| parihaka | psnr | 65 | 4.6970 | palette |  |  | 5.8155 |
| parihaka | psnr | 75 | 4.6970 | palette |  |  | 7.4792 |
| parihaka | psnr | 85 | 4.6970 | palette |  |  | 9.1704 |
| penobscot | linf | 0.00255 | 4.6735 | palette | yes | 7.3487 | 13.2689 |
| penobscot | linf | 0.0255 | 4.6735 | palette | yes | 7.3792 | 9.9470 |
| penobscot | linf | 0.255 | 4.6735 | palette | yes | 6.9236 | 6.6239 |
| penobscot | linf | 2.55 | 2.4939 | wavelet | yes | 3.5013 | 3.2869 |
| penobscot | psnr | 45 | 2.1055 | wavelet |  |  | 2.7578 |
| penobscot | psnr | 55 | 3.7580 | wavelet |  |  | 4.3471 |
| penobscot | psnr | 65 | 4.6735 | palette |  |  | 6.1345 |
| penobscot | psnr | 75 | 4.6735 | palette |  |  | 7.7290 |
| penobscot | psnr | 85 | 4.6735 | palette |  |  | 9.4876 |
| salt | linf | 0.02982 | 0.9545 | palette | yes | 2.0989 | 6.8530 |
| salt | linf | 0.2982 | 0.9545 | palette | yes | 1.7457 | 4.0716 |
| salt | linf | 2.982 | 0.3442 | palette | yes | 0.4026 | 1.9607 |
| salt | linf | 29.82 | 0.0987 | palette | yes | 0.2421 | 0.8230 |
| salt | psnr | 45 | 0.0964 | palette |  |  | 0.6864 |
| salt | psnr | 55 | 0.1925 | palette |  |  | 1.1718 |
| salt | psnr | 65 | 0.3442 | palette |  |  | 1.7536 |
| salt | psnr | 75 | 0.6764 | palette |  |  | 2.4419 |
| salt | psnr | 85 | 0.9545 | palette |  |  | 3.5739 |

### BD-rate @PSNR (conservative on lossless palette points)

| volume | SPERR | SZ3 | ZFP | TTHRESH | MGARD |
|---|---:|---:|---:|---:|---:|
| f3 | −10.3% | −20.2% | −22.5% | **+14.4%** | −32.7% |
| overthrust | −61.3% | −37.6% | −62.1% | −81.8% | −66.4% |
| parihaka | −33.6% | −30.6% | −53.1% | −39.7% | −50.9% |
| penobscot | −33.2% | −34.6% | −53.2% | −36.5% | −48.0% |
| salt | −74.8% | −10.6% | −72.3% | −92.0% | −56.3% |

### BD-rate @L-inf (same tau)

| volume | SPERR | SZ3 | ZFP | MGARD | tnvol bound ok | CPU vs SPERR |
|---|---:|---:|---:|---:|---|---:|
| f3 | −12.8% | −13.9% | −57.0% | −47.8% | yes | 1.64× |
| overthrust | −58.9% | −20.9% | −81.0% | −69.2% | yes | 0.86× |
| parihaka | −42.7% | −27.6% | −67.6% | −54.0% | yes | 1.52× |
| penobscot | −44.5% | −34.0% | −67.7% | −57.4% | yes | 1.12× |
| salt | −82.6% | −41.6% | −91.8% | −76.4% | yes | 0.85× |

### F3 vs TTHRESH at low rate (point, not BD-rate)

TTHRESH's official "PSNR 45" run on F3 achieved 52.91 dB at 1.843 bpv, so it
is not a matched-quality point. Interpolating TTHRESH's measured curve to
PSNR 45.0 dB gives about 0.951 bpv. r10's 0.7878 bpv at requested PSNR 45
(measured 45.27 dB on a crop check) is then **17.2% smaller than TTHRESH at
the same quality**, and 42.8% smaller than SPERR at the same requested point
(1.3773 bpv). This does not overturn the F3 BD-rate of +14.4% vs TTHRESH,
because r10's lead shrinks at high quality.

## Losses (summary)

- F3 vs TTHRESH: +14.4% BD-rate @PSNR.
- Sparse cloud / hydrometeor fields vs SZ3: `h_hurr_CLOUDlog` (+67.1% / +70.3%),
  `v_hurr_QICElog` (+63.6% / +73.8%), `w_hurr_QCLOUDlog` (+76.3% / +92.2%),
  plus smaller losses on `v_hurr_QRAINlog`, `v_hurr_QGRAUPlog`, and
  `h_hurr_QSNOWlog` at L-inf.
- `h_qmc` vs TTHRESH: +28.0% BD-rate @PSNR.
- Two tiny MGARD PSNR losses: `h_hurr_CLOUDlog` +0.7%, `v_hurr_QICElog` +0.7%.
