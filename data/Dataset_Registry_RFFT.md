# Dataset Registry for R-FFT-539.9

Created: 2026-07-27

## Binding Rule
Every series must be run free-first with phase-scramble nulls. discovery_claim_allowed only for free non-scrambled runs. Preconditioned is secondary labeled compatibility only.

## Staged Datasets (synthetic suite)
- pure_noise, pure_G4_tone, G4_plus_noise, off_target_600s, AM_G4, coloured_noise, multi_harmonic
- Sampling: 1 Hz, ~40 cycles of 539.9 s
- Formats: .npy (dict with t/signal/fs) and .csv
- Location (sandbox / local mirror): datasets/synthetic/

## Internal demo
- multi_k_seed_residuals: short vector from pure-even multi-k A5 seed coefficients (Category A arithmetic demo only)

## Public datasets (loader-only – fetch on local machine)
- LIGO O3 second-trend auxiliary channels (GWpy / NDS2, nds.gwosc.org)
- GWOSC calibrated strain segments (downsample to <=1 Hz)
- NANOGrav 15-yr / 12.5-yr residual series (Zenodo / nanograv.org)
- EHT M87* multi-epoch and DESI voids deferred (sparse or spatial)

## Loaders
See datasets/loaders/public_loaders.py and run_synthetics_free_first.py for stubs.

Free-first discipline is mandatory for any discovery claim.
