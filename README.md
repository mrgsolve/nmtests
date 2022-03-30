# nmtests

Benchmark tests against nonmem. These aren't tests of speed or performance, but
rather tests comparing different behaviors. The primary interest is things like
dosing events, steady state, handling time-varying covariates. 

The most-complete test document is in 
[dosing/dosing-vignette.md](dosing/dosing-vignette.md).

## Results

### Dosing vignette

- `dosing/results/1001.csv` - mrgsolve / nonmem comparison
- `dosing/results/1001-sims.csv` - mrgsolve simulated data
- `dosing/results/1001-nm.csv` - nonmem simulated data
- `dosing/results/1001.json` - latest run meta data
- `dosing/results/1001.pdf` - plots comparing nm and mrgsolve

## Requirements

- pkgr
- PsN
- NONMEM

## Installation

- Install [pkgr](https://github.com/metrumresearchgroup/pkgr)
- Install packages
  - `bash$ pkgr install`

## Run

All R vignette scripts run from the directory where the reside. 

- Run `dosing/dosing-vignette.R`
  - On the command line: `bash$ Make nmtest`
  - Via R studio: knit or run `dosing/dosing-vignette.R`
