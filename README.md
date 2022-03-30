# nmtests

Benchmark tests against nonmem. These aren't tests of speed or performance, but
rather tests comparing different behaviors. The primary interest is things like
dosing events, steady state, handling time-varying covariates. 

The most-complete test document is in 
[dosing/dosing-vignette.md](dosing/dosing-vignette.md). This is usually what
users are primarily interested in checking. 

There is another vignette testing time-varying covariates. The output is 
located in 
[time-varying/time-varying-vignette.md](time-varying/time-varying-vignette.md).


## Results

### Dosing vignette

- [dosing/dosing-vignette.md](dosing/dosing-vignette.md) - the completed vignette
- `dosing/results/1001.csv` - mrgsolve / nonmem comparison
- `dosing/results/1001-sims.csv` - mrgsolve simulated data
- `dosing/results/1001-nm.csv` - nonmem simulated data
- `dosing/results/1001.json` - latest run meta data
- `dosing/results/1001.pdf` - plots comparing nm and mrgsolve

### Time-varying covariate vignette

- [time-varying/time-varying-vignette.md](time-varying/time-varying-vignette.md) - the completed vignette
- `time-varying/results/2001.csv` - mrgsolve / nonmem comparison
- `time-varying/results/2001-sims.csv` - mrgsolve simulated data
- `time-varying/results/2001-nm.csv` - nonmem simulated data
- `time-varying/results/2001.json` - latest run meta data

## Requirements

- pkgr
  - Required R packages can be found in the [pkgr.yml](pkgr.yml) file
- PsN
- NONMEM

## Package installation

- Install [pkgr](https://github.com/metrumresearchgroup/pkgr)
- Install packages
  - `bash$ pkgr install`
- Any issues installing or using pkgr should be resolved with the pkgr 
  maintainers [here](https://github.com/metrumresearchgroup/pkgr/issues)

## Run

All R vignette scripts run from the directory where the reside. The vignette
simulates using both mrgsolve and NONMEM and requires NONMEM and PsN install; 
it is up to the user to install both and troubleshoot any issues with those
tools. There is also a `check` script that simulates with mrgsolve and compares
against previously simulated NONMEM output.

- Run `dosing/dosing-vignette.R`
  - On the command line: `bash$ make dosing`
  - To check: `bash$ make check-dosing`
  - Via R studio: knit or run `dosing/dosing-vignette.R`

- Run `time-varying/time-varying-vignette.R`
  - On the command line: `bash$ make tv`
  - To check: `bash$ make check-tv`
  - Via R studio: knit or run `time-varying/time-varying-vignette.R`
