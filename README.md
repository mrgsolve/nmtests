# nmtests

Benchmark tests against nonmem.


The most-complete test document is [here](nmtest9.md).

## Results
- `results/1001.csv` - mrgsolve / nonmem comparison
- `results/1001-sims.csv` - mrgsolve simulated data
- `results/1001-nm.csv` - nonmem simulated data
- `results/1001.json` - latest run meta data
- `results/1001.pdf` - plots comparing nm and mrgsolve

## Requirements

- pkgr
- PsN
- NONMEM

## Installation

- Install [pkgr](https://github.com/metrumresearchgroup/pkgr)
- Install packages
  - `bash$ pkgr install`

## Run
- Run `nmtest9.R`
  - On the command line: `bash$ Make nmtest`
  - Via R studio: knit or run `nmtest9.R`
