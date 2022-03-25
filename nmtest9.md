Bemchmark test with mrgsolve and NONMEM
================
Metrum Research Group

-   [mrgsolve package version](#mrgsolve-package-version)
-   [Introduction](#introduction)
-   [Setup](#setup)
-   [Functions](#functions)
    -   [Save `mrgsim` output as a `nonmem` input data
        set](#save-mrgsim-output-as-a-nonmem-input-data-set)
    -   [Save the `nonmem` input data
        set](#save-the-nonmem-input-data-set)
    -   [Run `nonmem`](#run-nonmem)
    -   [Read in `nonmem` simulation
        results](#read-in-nonmem-simulation-results)
    -   [Simulate a scenario with
        `mrsim()`](#simulate-a-scenario-with-mrsim)
-   [The `mrgsim()` model](#the-mrgsim-model)
-   [Assemble the scenarios](#assemble-the-scenarios)
-   [Simulate with `nonmem`](#simulate-with-nonmem)
-   [Numeric Summary](#numeric-summary)
    -   [Overall](#overall)
    -   [Summary by scenario number](#summary-by-scenario-number)
-   [Results](#results)
    -   [1: Bolus with additional](#1-bolus-with-additional)
    -   [2: Bolus with lag time and
        bioav](#2-bolus-with-lag-time-and-bioav)
    -   [3: Infusion with additional](#3-infusion-with-additional)
    -   [4: Infusion with bioav factor](#4-infusion-with-bioav-factor)
    -   [5: Infusion with bioav factor and
        dur](#5-infusion-with-bioav-factor-and-dur)
    -   [6: Infusion doses to depot, with
        additional](#6-infusion-doses-to-depot-with-additional)
    -   [7: Infusion doses, with additional and lag
        time](#7-infusion-doses-with-additional-and-lag-time)
    -   [8: Infusion doses, with lag time and bioav
        factor](#8-infusion-doses-with-lag-time-and-bioav-factor)
    -   [9: Infusion doses, with lag time and bioav
        factor](#9-infusion-doses-with-lag-time-and-bioav-factor)
    -   [10: Infusion doses at steady-state, with lag time and bioav
        factor](#10-infusion-doses-at-steady-state-with-lag-time-and-bioav-factor)
    -   [11: Infusion doses, with lag time and bioav
        factor](#11-infusion-doses-with-lag-time-and-bioav-factor)
    -   [12: Infusion doses at steady state, II &lt; DUR, no bioav
        factor](#12-infusion-doses-at-steady-state-ii--dur-no-bioav-factor)
    -   [13: Infusion doses at steady state where II == DUR, with bioav
        factor](#13-infusion-doses-at-steady-state-where-ii--dur-with-bioav-factor)
    -   [14: Infusion doses at steady state, where II ==
        DUR](#14-infusion-doses-at-steady-state-where-ii--dur)
    -   [15: Bolus doses at steady state, with bioav factor and lag
        time](#15-bolus-doses-at-steady-state-with-bioav-factor-and-lag-time)
    -   [16: Bolus doses with lag time and bioavability
        factor](#16-bolus-doses-with-lag-time-and-bioavability-factor)
    -   [17: Bolus then infusion](#17-bolus-then-infusion)
    -   [18: Infusion with modeled duration, lag time, and bioav
        factor](#18-infusion-with-modeled-duration-lag-time-and-bioav-factor)
    -   [19: Infusion with modeled duration, at steady state with bioav
        factor](#19-infusion-with-modeled-duration-at-steady-state-with-bioav-factor)
    -   [20: Reset and dose (EVID 4) with
        additional](#20-reset-and-dose-evid-4-with-additional)
    -   [21: Reset (EVID 3) with
        additional](#21-reset-evid-3-with-additional)
-   [Control stream](#control-stream)
-   [Session Info](#session-info)

# mrgsolve package version

``` r
packageVersion("mrgsolve")
```

    ## [1] '1.0.3'

# Introduction

This document runs simulations from a pharmacokinetic model using both
NONMEM and mrgsolve and compares the results.

All of the relevant code is presented so that the user can trace how the
simulations are performed. The complete source code can be viewed
[here](nmtest9.R).

The bottom line results are presented in graphical form [here](#results)
and numeric form [here](#numeric-summary).

# Setup

``` r
Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")
```

``` r
mrgsolve.loc <- NULL
library(mrgsolve, lib.loc = mrgsolve.loc)
```

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
library(ggplot2)
library(purrr)
```

    ## 
    ## Attaching package: 'purrr'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     transpose

``` r
library(tidyr)
library(jsonlite)
```

    ## 
    ## Attaching package: 'jsonlite'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     flatten

``` r
library(tools)
library(parallel)
```

``` r
stopifnot(grepl("PsN", system("execute --version", intern=TRUE)))
```

# Functions

These functions assemble data sets, run simulations, and gather outputs.
All scenarios are handled in exactly the same way.

## Save `mrgsim` output as a `nonmem` input data set

``` r
to_data_set <- function(x, id = NULL) {
  x <- as.data.frame(x)
  x <- mutate(x, C = '.', DV = '.', cmt = if_else(cmt==0, 2, cmt))
  x <- select(x, "C", everything())
  if(is.numeric(id)) x <- mutate(x,ID = id)
  x
}
```

## Save the `nonmem` input data set

``` r
sv <- function(x,file) {
  fwrite(file = file, quote = FALSE, x, na = '.')
}
```

## Run `nonmem`

``` r
run <- function(number) {
  #execute -model_subdir -silent -directory=run1001 model/1001.ctl
  directory <- paste0("-directory=", number)
  target <- paste0("model/", number, ".ctl")
  args <- c("-silent", "-model_subdir", directory, target)
  system2("execute", args=args)
  return(tabread(number))
}
```

## Read in `nonmem` simulation results

``` r
tabread <- function(number) {
  tab <- file.path(number, "TAB")
  if(file.exists(tab)) return(fread(tab, skip=1))
  stop("the run failed")
}
```

## Simulate a scenario with `mrsim()`

``` r
sim <- function(x, e,...) {
  carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")
  mrgsim(x, events = e, carry_out = carry, digits = 5, recsort = 3, ...) 
}
```

# The `mrgsim()` model

``` r
mod <- mread_cache("model/1001.mod")
```

    . Building 1001_mod ... done.

``` r
mod <- update(mod, end = 130, delta = 1)
```

# Assemble the scenarios

There is a lot of code here. See the [results](#results) section to see
input data objects next to simulated data output from mrgsolve and
NONMEM.

-   Doses into `cmt` 2 are intravascular and doses into `cmt` 1 are
    extravascular
-   `LAGT` sets the dosing lag time
-   `BIOAV` sets the bioavailability fraction

``` r
env <- new.env()
env$ev <- list()
env$descr <- list()
push_back <- function(env, ev, descr) {
  n <- length(env$ev)+1
  m <- length(env$descr)+1
  env$ev[[n]] <- ev
  env$descr[[m]] <- descr
}
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3) 
push_back(env,ev, "Bolus with additional")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, LAGT = 12.13, BIOAV = 2.23, cmt = 2) 
push_back(env, ev,"Bolus with lag time and bioav")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, cmt = 2) 
push_back(env,ev,"Infusion with additional")
```

``` r
ev <- ev(amt = 480, ii = 0, addl = 0, rate = 10, cmt = 2, BIOAV = 0.5) 
push_back(env,ev,"Infusion with bioav factor")
```

``` r
ev <- ev(amt = 480, ii = 0, addl = 0, rate = -2, cmt = 2, BIOAV = 0.5, MODE = 2, DUR2 = 48) 
push_back(env,ev,"Infusion with bioav factor and dur")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/12, cmt = 1) 
push_back(env,ev,"Infusion doses to depot, with additional")
```

``` r
ev <- ev(amt = 100, ii = 24, addl=3, rate = 100/10, LAGT = 4.15, cmt = 2) 
push_back(env,ev,"Infusion doses, with additional and lag time")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.25, BIOAV = 0.412, cmt = 2) 
push_back(env,ev,"Infusion doses, with lag time and bioav factor")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.16, BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses, with lag time and bioav factor")
```

``` r
ev <- ev(amt = 100, ii = 12, addl = 4, rate = 100/50, BIOAV = 0.812, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady-state, with lag time and bioav factor")
```

``` r
ev <- ev(amt = 100, ii = 12, addl = 3, rate = 100/50, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses, with lag time and bioav factor")
```

``` r
ev <- ev(amt = 100, ii = 6, addl = 12, rate = signif(100/12,5), ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady state, II < DUR, no bioav factor")
```

``` r
ev <- ev(amt = 100, ii = 10, addl = 8, rate = 0.412*100/10,  BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady state where II == DUR, with bioav factor")
```

``` r
ev <- ev(amt = 100, ii = 10, addl = 8, rate = 100/10, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady state, where II == DUR")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3,  LAGT = 4, BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Bolus doses at steady state, with bioav factor and lag time")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3,  LAGT = 5, BIOAV = 0.412, cmt = 2) 
push_back(env,ev,"Bolus doses with lag time and bioavability factor")
```

``` r
ev <- 
  ev(amt = 100, cmt = 2, LAGT = 1) + 
  ev(time = 13, amt = 50, ii = 24, addl = 2, rate = 24)
push_back(env,ev,"Bolus then infusion")
```

``` r
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, LAGT = 5, BIOAV = 0.61)
push_back(env,ev,"Infusion with modeled duration, lag time, and bioav factor")
```

``` r
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, ss = 1, BIOAV = 0.61)
push_back(env,ev,"Infusion with modeled duration, at steady state with bioav factor")
```

``` r
ev <- 
  ev(amt = 100, ii = 12, addl = 2, rate = 50, BIOAV = 0.61) + 
  ev(amt = 120, evid = 4, time = 50, BIOAV = 0.5, ii = 12, addl = 3)
push_back(env,ev,"Reset and dose (EVID 4) with additional")
```

``` r
ev <- 
  ev(amt = 100, ii = 12, addl = 3, rate = 50, BIOAV = 0.61) + 
  ev(amt = 0, evid = 3, time = 50, cmt = 2, BIOAV=1) + 
  ev(amt = 120, ii = 16, addl = 2, time = 54, BIOAV=1)
push_back(env,ev,"Reset (EVID 3) with additional")
```

``` r
ev <- 
  ev(amt = 100, ii = 24, addl = 3, ss = 1)  + 
  ev(amt = 50,  ii = 24, addl = 3, ss = 2, time = 12)
push_back(env,ev,"Steady state 1 and 2")
```

``` r
ev <- ev(amt = 0, rate = 100,  ss=1)
push_back(env,ev,"Steady state infusion")
```

``` r
update_id <- function(ev,id) mutate(ev, ID = id)
```

``` r
runs <- tibble(ev = env$ev, descr = env$descr)
runs <- mutate(runs, ID = seq(n()))
runs <- mutate(runs, ev = map2(ev, ID, update_id))
runs <- mutate(runs, sims = mclapply(ev, sim, x = mod))
runs <- mutate(runs, data = map(sims, to_data_set))
```

``` r
data <- bind_rows(runs$data)
data$CP <- NULL

sv(data, "data/1001.csv")
```

# Simulate with `nonmem`

``` r
nm <- run(1001)
```

# Numeric Summary

Look at the difference between simulated values from mrgsolve and
NONMEM.

``` r
runs <- mutate(runs, nm = split(nm, nm$ID))

runs <- mutate(
  runs, 
  comp = map2(nm, sims, .f = function(nonmem, mrgsolve) {
    tibble(
      ID = nonmem$ID, 
      time = mrgsolve$time, 
      MRGSIM = mrgsolve$CP, 
      NONMEM = nonmem$CP, 
      diff = MRGSIM-NONMEM
    )  
  })
)

comp <- select(runs, comp) %>% unnest(cols = c(comp))
sims <- pull(runs, sims) %>% lapply(as_tibble) %>% rbindlist()
fwrite(x = comp, file = "results/1001.csv")
fwrite(x = sims, file = "results/1001-sims.csv")
fwrite(x = nm, file = "results/1001-nm.csv")

meta <- list()
meta$md5 <- list()
meta$md5$data <- md5sum("data/1001.csv")
meta$md5$ctl <- md5sum("model/1001.ctl")
meta$md5$mod <- md5sum("model/1001.mod")
meta$md5$result <- md5sum("results/1001.csv")
meta$md5$sims <- md5sum("results/1001-sims.csv")
meta$md5$sims <- md5sum("results/1001-nm.csv")
meta$date <- date()
meta$user <- Sys.info()[["user"]]
write_json(
  x = meta, 
  path = "results/1001.json",
  pretty = TRUE
)
```

## Overall

This is the `nonmem` minus `mrgsim()` summary

``` r
summary(comp$diff)
```

    .    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    .       0       0       0       0       0       0

## Summary by scenario number

`diff` is the simulated `CP` from `nonmem` minus the simulated `CP` from
`mrgsim()`

``` r
group_by(comp, ID) %>% 
  summarise(mean = mean(diff), max = max(diff), min = min(diff)) %>% 
  as.data.frame()
```

    .    ID mean max min
    . 1   1    0   0   0
    . 2   2    0   0   0
    . 3   3    0   0   0
    . 4   4    0   0   0
    . 5   5    0   0   0
    . 6   6    0   0   0
    . 7   7    0   0   0
    . 8   8    0   0   0
    . 9   9    0   0   0
    . 10 10    0   0   0
    . 11 11    0   0   0
    . 12 12    0   0   0
    . 13 13    0   0   0
    . 14 14    0   0   0
    . 15 15    0   0   0
    . 16 16    0   0   0
    . 17 17    0   0   0
    . 18 18    0   0   0
    . 19 19    0   0   0
    . 20 20    0   0   0
    . 21 21    0   0   0
    . 22 22    0   0   0
    . 23 23    0   0   0

``` r
comp_plot <- function(comp) {
  id <- comp$ID[1]
  title <- paste0("ID: ", id, "; line: mrgsolve, point: NONMEM")
  ggplot(data = comp) + 
    geom_point(aes(time, NONMEM),color = "firebrick") + 
    geom_line(aes(time, MRGSIM, group = ID)) +
    ylab("Simulated value") + xlab("Time") + 
    ggtitle(label = NULL, subtitle = title) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(0, 130, 24))  
}

runs <- mutate(runs, plot = map(comp, comp_plot))
```

# Results

## 1: Bolus with additional

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,7]> <gg>

## 2: Bolus with lag time and bioav

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 3: Infusion with additional

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,8]> <gg>

## 4: Infusion with bioav factor

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 5: Infusion with bioav factor and dur

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,11]> <gg>

## 6: Infusion doses to depot, with additional

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,8]> <gg>

## 7: Infusion doses, with additional and lag time

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 8: Infusion doses, with lag time and bioav factor

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,10]> <gg>

## 9: Infusion doses, with lag time and bioav factor

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,11]> <gg>

## 10: Infusion doses at steady-state, with lag time and bioav factor

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,10]> <gg>

## 11: Infusion doses, with lag time and bioav factor

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 12: Infusion doses at steady state, II &lt; DUR, no bioav factor

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 13: Infusion doses at steady state where II == DUR, with bioav factor

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,10]> <gg>

## 14: Infusion doses at steady state, where II == DUR

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 15: Bolus doses at steady state, with bioav factor and lag time

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,10]> <gg>

## 16: Bolus doses with lag time and bioavability factor

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 17: Bolus then infusion

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 18: Infusion with modeled duration, lag time, and bioav factor

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,12]> <gg>

## 19: Infusion with modeled duration, at steady state with bioav factor

    . # A tibble: 1 × 2
    .   ev        plot        
    .   <list>    <named list>
    . 1 <ev[,12]> <gg>

## 20: Reset and dose (EVID 4) with additional

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

## 21: Reset (EVID 3) with additional

    . # A tibble: 1 × 2
    .   ev       plot        
    .   <list>   <named list>
    . 1 <ev[,9]> <gg>

# Control stream

``` r
writeLines(readLines("model/1001/1001.lst"))
```

       Thu Mar 24 11:19:21 EDT 2022
       $PROB RUN# 101
       
       $INPUT C ID TIME EVID AMT CMT SS II ADDL RATE LAGT MODE DUR2 RAT2 BIOAV DV
       
       $DATA ../../data/1001.csv IGNORE=C
       
       $SUBROUTINES ADVAN2 TRANS2
       
       $PK
       
       TVCL=THETA(1)
       CL=TVCL*EXP(ETA(1))
       
       TVV2=THETA(2)
       V=TVV2*EXP(ETA(2))
       
       TVKA=THETA(3)
       KA=TVKA*EXP(ETA(3))
       
       ALAG2 = LAGT
       F2 = BIOAV
       
       IF(MODE.EQ.1) R2 = RAT2
       IF(MODE.EQ.2) D2 = DUR2
       
       $ERROR
       IPRED=A(2)/(V/1000)
       Y=IPRED*EXP(ERR(1))
       
       CP = IPRED
       
       $THETA
       (1.1,   FIX) ;; CL
       (20,  FIX) ;; V
       (1.5, FIX) ;; KA
       
       $OMEGA
       0.0 FIX
       0.0 FIX
       0.0 FIX
       
       $SIGMA
       0.00 FIX
       
       $TABLE FILE=TAB ID TIME CP NOPRINT ONEHEADER NOAPPEND
       
       $SIMULATION (2674474) ONLYSIMULATION
       
       
       NM-TRAN MESSAGES
         
        WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
                    
        (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
       
       License Registered to: Metrum Research Group (with RADAR5NM)
       Expiration Date:    14 JUL 2022
       Current Date:       24 MAR 2022
       Days until program expires : 110
       1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.4.4
        ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
        CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
        AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
        PERFORMED BY NOUS INFOSYSTEMS.
       
        PROBLEM NO.:         1
        RUN# 101
       0DATA CHECKOUT RUN:              NO
        DATA SET LOCATED ON UNIT NO.:    2
        THIS UNIT TO BE REWOUND:        NO
        NO. OF DATA RECS IN DATA SET:     3041
        NO. OF DATA ITEMS IN DATA SET:  17
        ID DATA ITEM IS DATA ITEM NO.:   2
        DEP VARIABLE IS DATA ITEM NO.:  16
        MDV DATA ITEM IS DATA ITEM NO.: 17
       0INDICES PASSED TO SUBROUTINE PRED:
          4   3   5  10   7   8   6   0   0   0   9
       0LABELS FOR DATA ITEMS:
        C ID TIME EVID AMT CMT SS II ADDL RATE LAGT MODE DUR2 RAT2 BIOAV DV MDV
       0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
        CP
       0FORMAT FOR DATA:
        (E2.0,E3.0,E4.0,E2.0,E4.0,2E2.0,2E3.0,E17.0,E6.0,E2.0,4E6.0,1F2.0)
       
        TOT. NO. OF OBS RECS:     3013
        TOT. NO. OF INDIVIDUALS:       23
       0LENGTH OF THETA:   3
       0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
       0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
       0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
       0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
       0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
       0INITIAL ESTIMATE OF THETA:
        LOWER BOUND    INITIAL EST    UPPER BOUND
         0.1100E+01     0.1100E+01     0.1100E+01
         0.2000E+02     0.2000E+02     0.2000E+02
         0.1500E+01     0.1500E+01     0.1500E+01
       0INITIAL ESTIMATE OF OMEGA:
        0.0000E+00
        0.0000E+00   0.0000E+00
        0.0000E+00   0.0000E+00   0.0000E+00
       0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
       0INITIAL ESTIMATE OF SIGMA:
        0.0000E+00
       0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
       0SIMULATION STEP OMITTED:    NO
        OBJ FUNC EVALUATED:         NO
        ORIGINAL DATA USED ON EACH NEW SIMULATION:         NO
        SEEDS RESET ON EACH NEW SUPERSET ITERATION:        YES
       0SIMULATION RANDOM METHOD SELECTED (RANMETHOD): 4U
       SEED   1 RESET TO INITIAL: YES
        SOURCE   1:
          SEED1:       2674474   SEED2:             0   PSEUDO-NORMAL
       0WARNING: NO. OF OBS RECS IN INDIVIDUAL REC NO.      1 (IN INDIVIDUAL
        REC ORDERING) EXCEEDS ONE WHILE INITIAL ESTIMATE OF WITHIN INDIVIDUAL VARIANCE IS ZERO
       0WARNING: NO. OF OBS RECS IN INDIVIDUAL REC NO.      2 (IN INDIVIDUAL
        REC ORDERING) EXCEEDS ONE WHILE INITIAL ESTIMATE OF WITHIN INDIVIDUAL VARIANCE IS ZERO
       0WARNING: NO. OF OBS RECS IN INDIVIDUAL REC NO.      3 (IN INDIVIDUAL
        REC ORDERING) EXCEEDS ONE WHILE INITIAL ESTIMATE OF WITHIN INDIVIDUAL VARIANCE IS ZERO
       0WARNING: NO. OF OBS RECS IN INDIVIDUAL REC NO.      4 (IN INDIVIDUAL
        REC ORDERING) EXCEEDS ONE WHILE INITIAL ESTIMATE OF WITHIN INDIVIDUAL VARIANCE IS ZERO
       0WARNING: NO. OF OBS RECS IN INDIVIDUAL REC NO.      5 (IN INDIVIDUAL
        REC ORDERING) EXCEEDS ONE WHILE INITIAL ESTIMATE OF WITHIN INDIVIDUAL VARIANCE IS ZERO
       0TABLES STEP OMITTED:    NO
        NO. OF TABLES:           1
        SEED NUMBER (SEED):    11456
        RANMETHOD:             3U
        MC SAMPLES (ESAMPLE):    300
        WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE
       0-- TABLE   1 --
       0RECORDS ONLY:    ALL
       04 COLUMNS APPENDED:    NO
        PRINTED:                NO
        HEADERS:               ONE
        FILE TO BE FORWARDED:   NO
        FORMAT:                S1PE11.4
        LFORMAT:
        RFORMAT:
        FIXED_EFFECT_ETAS:
       0USER-CHOSEN ITEMS:
        ID TIME CP
       1DOUBLE PRECISION PREDPP VERSION 7.4.4
       
        ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
       0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
       0BASIC PK PARAMETERS (AFTER TRANSLATION):
          ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
          ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3
       
        TRANSLATOR WILL CONVERT PARAMETERS
        CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
       0COMPARTMENT ATTRIBUTES
        COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                                STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
           1         DEPOT        OFF        YES        YES        YES        NO
           2         CENTRAL      ON         NO         YES        NO         YES
           3         OUTPUT       OFF        YES        NO         NO         NO
       1
        ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
        COMPT. NO.                             INDICES
                     SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                                FRACTION    RATE        DURATION    LAG
           1            *           *           *           *           *
           2            *           5           6           7           4
           3            *           -           -           -           -
                    - PARAMETER IS NOT ALLOWED FOR THIS MODEL
                    * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
                      WILL DEFAULT TO ONE IF APPLICABLE
       0DATA ITEM INDICES USED BY PRED ARE:
          EVENT ID DATA ITEM IS DATA ITEM NO.:      4
          TIME DATA ITEM IS DATA ITEM NO.:          3
          DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   5
          DOSE RATE DATA ITEM IS DATA ITEM NO.:    10
          STEADY STATE DATA ITEM IS DATA ITEM NO.:  7
          INTERVAL DATA ITEM IS DATA ITEM NO.:      8
          ADDL. DOSES DATA ITEM IS DATA ITEM NO.:   9
          COMPT. NO. DATA ITEM IS DATA ITEM NO.:    6
       
       0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
        PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
       0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
       0ERROR SUBROUTINE INDICATES THAT DERIVATIVES OF COMPARTMENT AMOUNTS ARE USED.
       1
        SIMULATION STEP PERFORMED
        SOURCE  1:
           SEED1:     918246936   SEED2:             0
        Elapsed simulation  time in seconds:     0.00
        ESTIMATION STEP OMITTED:                 YES
        Elapsed finaloutput time in seconds:     0.04
        #CPUT: Total CPU Time in Seconds,        0.060
       Stop Time:
       Thu Mar 24 11:19:27 EDT 2022

# Session Info

``` r
options(width = 120)
sessionInfo()
```

    . R version 4.1.1 (2021-08-10)
    . Platform: x86_64-pc-linux-gnu (64-bit)
    . Running under: Ubuntu 18.04.5 LTS
    . 
    . Matrix products: default
    . BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    . LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
    . 
    . locale:
    .  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    .  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    .  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    . 
    . attached base packages:
    . [1] parallel  tools     stats     graphics  grDevices datasets  utils     methods   base     
    . 
    . other attached packages:
    . [1] jsonlite_1.8.0    tidyr_1.2.0       purrr_0.3.4       ggplot2_3.3.5     data.table_1.14.2 dplyr_1.0.8      
    . [7] mrgsolve_1.0.3   
    . 
    . loaded via a namespace (and not attached):
    .  [1] Rcpp_1.0.8.3     pillar_1.7.0     compiler_4.1.1   highr_0.9        digest_0.6.29    evaluate_0.15   
    .  [7] lifecycle_1.0.1  tibble_3.1.6     gtable_0.3.0     pkgconfig_2.0.3  rlang_1.0.2      cli_3.2.0       
    . [13] yaml_2.3.5       xfun_0.30        fastmap_1.1.0    withr_2.5.0      stringr_1.4.0    knitr_1.37      
    . [19] generics_0.1.2   vctrs_0.3.8      grid_4.1.1       tidyselect_1.1.2 glue_1.6.2       R6_2.5.1        
    . [25] fansi_1.0.2      rmarkdown_2.13   magrittr_2.0.2   scales_1.1.1     ellipsis_0.3.2   htmltools_0.5.2 
    . [31] colorspace_2.0-3 renv_0.14.0      utf8_1.2.2       stringi_1.7.6    munsell_0.5.0    crayon_1.5.0
