# Dataset Registry for R-FFT-539.9 (Expanded Observatories Edition)

Created / updated: 2026-07-27

## Binding Rule (mandatory)

Every series must be run free-first with phase-scramble nulls.
discovery_claim_allowed = True only for free (precondition=False), non-scrambled runs.
Preconditioned lock is secondary and must be labelled compatibility-only.
The refined_period under precondition=True is forced toward the assumed G4 by construction; it is not an independent discovery statistic.

## Staged Datasets (already written)

### Synthetic suite (1 Hz, ~40 cycles of 539.9 s)
- pure_noise, pure_G4_tone, G4_plus_noise, off_target_600s, AM_G4, coloured_noise, multi_harmonic
- Paths: datasets/synthetic/*.npy and *.csv
- free_first_required: True

### Internal quantum-algebraic placeholder
- multi_k_seed_residuals (short vector from pure-even multi-k A5 catalogue)
- Path: datasets/internal/
- free_first_required: True
- Note: Category A arithmetic demo only; not a long physical time series.

### Illustrative PTA residual sample
- example_pta_residual_sample.csv (synthetic microsecond-scale residual for loader testing)

## Observatories Section

### 1. LIGO (Hanford H1, Livingston L1)
- Best public access: GWOSC (https://gwosc.org) + GWpy / NDS2 (nds.gwosc.org)
- Preferred for ~540 s search: O3/O4 second-trend auxiliary channels (seismometers, accelerometers, RMS, wind) at ~1 s sampling. Multi-month continuous series exist.
- Strain: high-rate (4–16 kHz); must be heavily downsampled or restricted to quiet continuous segments.
- Free-first required: Yes. Always run free estimator + nulls first.
- Loader note: use TimeSeries.fetch(channel, start, end, host="nds.gwosc.org")

### 2. Virgo (V1)
- Best public access: same GWOSC / GWpy / NDS2 path as LIGO.
- Channels: V1:PEM-* or equivalent second-trend / environmental monitors where released.
- Sampling and free-first rules identical to LIGO.
- Loader note: substitute V1 channel names in the GWpy fetch example.

### 3. KAGRA (K1)
- Best public access: GWOSC open data releases that include KAGRA.
- Availability of long second-trend auxiliary is more limited than LIGO/Virgo; check current GWOSC catalogue.
- Free-first required: Yes.
- Loader note: K1 channel names via GWOSC timeline / GWpy when present.

### 4. GEO 600
- Best public access: limited open segments via GWOSC.
- Useful mainly for short comparison or detector-characterisation studies.
- Free-first required: Yes.

### 5. NANOGrav / PTA (and related IPTA, EPTA, PPTA data)
- Best public access: nanograv.org/science/data and Zenodo (15-year and 12.5-year residual / TOA releases).
- Data type: post-fit timing residuals and TOAs. Native sensitivity is nanohertz (periods of years).
- For R-FFT: convert residual tables to an evenly sampled residual vector (interpolate or bin carefully). Useful as long-baseline control domain or for super-harmonic checks; 539.9 s is far below the typical PTA Nyquist.
- Free-first required: Yes.
- Loader note: parse .tim / residual ASCII or HDF5, then np.interp or bin to regular grid.

### 6. Event Horizon Telescope (EHT) / M87*
- Best public access: CyVerse Data Commons, GitHub eventhorizontelescope repositories (visibility amplitudes, polarised data for 2017/2018/… campaigns), accompanying papers.
- Sampling: sparse VLBI epochs (hours to days between scans, years between campaigns). Not continuous high-cadence.
- Practical use for R-FFT: construct a proxy light-curve or polarisation time series from published monitoring points, then test free-first. Native data are better suited to multi-year polarity / morphology studies than a 9-minute continuous search.
- Free-first required: Yes.
- Loader note: extract time-stamped amplitude / polarisation values from public tables or UVFITS summaries, interpolate only if scientifically justified.

### 7. Other facilities (noted, deferred or limited)
- LISA: future; no long public continuous series yet.
- Planck / WMAP: CMB maps are spatial; time-ordered data exist but are specialised and not primary for 540 s searches.
- DESI / SDSS / BOSS: void and LSS catalogues are spatial. Stacked ISW profiles exist but are not long temporal sequences.
- Fermi / Swift / GRB monitors: light curves are event-triggered and usually short; selected long monitoring light curves may be usable after even sampling.
- Atomic-clock / precision timing networks: public residual series exist for some systems but must be located case-by-case.

## Free-first execution order (every observatory)

1. Load or construct the evenly sampled series.
2. Run free estimator (precondition=False) + phase-scramble / scrambled-G4 nulls.
3. Record free_T_hat, SNR, discovery_claim_allowed.
4. Only then optionally run preconditioned mode and label it compatibility-only.
5. Never present a preconditioned refined_period as unsupervised discovery of 539.9 s.

## Public bulk series status
Full multi-month LIGO/Virgo/KAGRA auxiliary archives, complete NANOGrav residual releases, and full EHT visibility sets remain loader-only (size limits). Use the concrete loader examples on a networked machine that already hosts scripts/r_fft_5399.py.
