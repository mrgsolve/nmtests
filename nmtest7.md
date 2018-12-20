Tests with NONMEM
================
Metrum Research Group, LLC

-   [Functions](#functions)
    -   [Save `mrgsim` output as a `nonmem` input data set](#save-mrgsim-output-as-a-nonmem-input-data-set)
    -   [Save the `nonmem` input data set](#save-the-nonmem-input-data-set)
    -   [Run `nonmem`](#run-nonmem)
    -   [Read in `nonmem` simulation results](#read-in-nonmem-simulation-results)
    -   [Simulate a scenario with `mrsim`](#simulate-a-scenario-with-mrsim)
-   [The `mrgsim` model](#the-mrgsim-model)
-   [Assemble the scenarios](#assemble-the-scenarios)
-   [Simulate with `nonmem`](#simulate-with-nonmem)
-   [Overall Summary](#overall-summary)
-   [Summary by RUN](#summary-by-run)
-   [Results](#results)
    -   [Bolus with additional](#bolus-with-additional)
    -   [Bolus with lag time and bioav](#bolus-with-lag-time-and-bioav)
    -   [Infusion with additional](#infusion-with-additional)
    -   [Infusion doses to depot, with additional](#infusion-doses-to-depot-with-additional)
    -   [Infusion doses, with additional and lag time](#infusion-doses-with-additional-and-lag-time)
    -   [Infusion doses, with lag time and bioav factor](#infusion-doses-with-lag-time-and-bioav-factor)
    -   [Infusion doses, with lag time and bioav factor](#infusion-doses-with-lag-time-and-bioav-factor-1)
    -   [Infusion doses at steady-state, with lag time and bioav factor](#infusion-doses-at-steady-state-with-lag-time-and-bioav-factor)
    -   [Infusion doses, with lag time and bioav factor](#infusion-doses-with-lag-time-and-bioav-factor-2)
    -   [Infusion doses at steady state, II &lt; DUR, no bioav factor](#infusion-doses-at-steady-state-ii-dur-no-bioav-factor)
    -   [Infusion doses at steady state where II == DUR, with bioav factor](#infusion-doses-at-steady-state-where-ii-dur-with-bioav-factor)
    -   [Infusion doses at steady state, where II == DUR](#infusion-doses-at-steady-state-where-ii-dur)
    -   [Bolus doses at steady state, with bioav factor and lag time](#bolus-doses-at-steady-state-with-bioav-factor-and-lag-time)
    -   [Bolus doses with lag time and bioavability factor](#bolus-doses-with-lag-time-and-bioavability-factor)
    -   [Bolus then infusion](#bolus-then-infusion)
    -   [Infusion with modeled duration, lag time, and bioav factor](#infusion-with-modeled-duration-lag-time-and-bioav-factor)
    -   [Infusion with modeled duration, at steady state with bioav factor](#infusion-with-modeled-duration-at-steady-state-with-bioav-factor)
    -   [Reset and dose (EVID 4) with additional](#reset-and-dose-evid-4-with-additional)
    -   [Reset (EVID 3) with additional](#reset-evid-3-with-additional)
    -   [Steady state 1 and 2](#steady-state-1-and-2)
-   [Control stream](#control-stream)
-   [Session Info](#session-info)

``` r
Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")
```

``` r
.libPaths("/data/Rlibs")
library(mrgsolve)
library(dplyr)
library(readr)
library(ggplot2)
library(parallel)
library(purrr)
library(tidyr)
```

``` r
carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")
```

Functions
=========

Save `mrgsim` output as a `nonmem` input data set
-------------------------------------------------

``` r
to_data_set <- function(x,id=NULL) {
  x <- as.data.frame(x)
  x <- mutate(x, C = '.', DV = '.', cmt = if_else(cmt==0, 2, cmt))
  x <- dplyr::select(x, "C", everything())
  if(is.numeric(id)) x <- mutate(x,ID = id)
  x
}
```

Save the `nonmem` input data set
--------------------------------

``` r
sv <- function(x,file) {
  write.csv(file = file, row.names = FALSE, quote = FALSE, x)
}
```

Run `nonmem`
------------

``` r
run <- function(number) {
  metrumrg::NONR(number, project = "model", 
                 command = "/opt/NONMEM/nm74/nmqual/autolog.pl", 
                 checkrunno=FALSE)
  return(tabread(number))
}
```

Read in `nonmem` simulation results
-----------------------------------

``` r
tabread <- function(number) {
  tab <- file.path("model", number, "TAB")
  if(file.exists(tab)) return(read_table(tab, skip=1))
  stop("the run failed")
}
```

Simulate a scenario with `mrsim`
--------------------------------

``` r
sim <- function(x, e,...) {
  mrgsim(x, events = e, carry.out = carry, digits = 5, ...) 
}

nmtest_plot <- function(x) {
  plot(x, CP~time)  
}
```

The `mrgsim` model
==================

``` r
code <- '
$SET req = ""
$PARAM CL = 1.1, V = 20, KA = 1.5
LAGT = 0, MODE = 0, DUR2 = 2, RAT2 = 10, BIOAV = 1

$PKMODEL cmt = "GUT CENT", depot = TRUE

$MAIN

F_CENT = BIOAV;
ALAG_CENT = LAGT;

if(MODE==1) R_CENT = RAT2;
if(MODE==2) D_CENT = DUR2;


$TABLE
capture DV = (CENT/(V/1000));
capture CP = DV;

$CAPTURE LAGT MODE DUR2 RAT2 BIOAV
'

mod <- mcode_cache("tests1", code)
```

    . Building tests1 ... done.

``` r
mod <- update(mod, end=130, delta = 1)
```

Assemble the scenarios
======================

``` r
env <- new.env()
env$ev <- list()
env$descr <- list()
push_back <- function(env,ev,descr) {
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
push_back(env, ev, "Bolus with lag time and bioav")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, cmt = 2) 
push_back(env,ev, "Infusion with additional")
```

``` r
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/12, cmt = 1) 
push_back(env,ev,"Infusion doses to depot, with additional")
```

``` r
ev4 <- ev(amt = 100, ii = 24, addl=3, rate = 100/10, LAGT = 4.15, cmt = 2) 
push_back(env,ev,"Infusion doses, with additional and lag time")
```

``` r
ev5 <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.25, BIOAV = 0.412, cmt = 2) 
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
ev <- ev(amt = 100, ii = 24, addl=3,  LAGT = 4, BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Bolus doses at steady state, with bioav factor and lag time")
```

``` r
ev <- ev(amt = 100, ii = 24, addl=3,  LAGT = 5, BIOAV = 0.412, cmt = 2) 
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
  ev(amt = 0, evid = 3, time = 50, cmt = 2) + 
  ev(amt = 120, ii = 16, addl = 2, time = 54)
push_back(env,ev,"Reset (EVID 3) with additional")
```

``` r
ev <- 
  ev(amt = 100, ii = 24, addl = 3, ss = 1)  + 
  ev(amt = 50,  ii = 24, addl = 3, ss = 2, time = 12)
push_back(env,ev,"Steady state 1 and 2")
```

``` r
update_id <- function(ev,id) mutate(ev, ID = id)

runs <- data_frame(ev = env$ev, descr = env$descr)
runs <- mutate(runs, ID = seq(n()))
runs <- mutate(runs,ev = map2(ev,ID, update_id))
runs <- mutate(runs, sims = mclapply(ev, sim, x = mod))

runs <- mutate(runs, data = map(sims, to_data_set))
```

``` r
data <- runs[["data"]] %>% bind_rows()

sv(data, "data/1001.csv")
```

Simulate with `nonmem`
======================

``` r
out <- run(1001)
```

    . Run 1001 complete.

    . NONR complete.

    . Parsed with column specification:
    . cols(
    .   ID = col_double(),
    .   TIME = col_double(),
    .   EVID = col_double(),
    .   CP = col_double(),
    .   IPRED = col_double(),
    .   PRED = col_double(),
    .   DV = col_double()
    . )

Overall Summary
===============

``` r
runs <- mutate(runs, out = split(out,out$ID))

runs <- mutate(
  runs, 
  comp = map2(out,sims, .f=function(out,sims) {
    data_frame(
      ID = out$ID, 
      time = sims$time, 
      MRGSIM = sims$CP, 
      NONMEM = out$CP, 
      diff = MRGSIM-NONMEM)  
  })
)

comp <- runs %>% select(ID,comp) %>% unnest()
```

This is the `nonmem` minus `mrgsim` summary

``` r
summary(comp$diff)
```

    .    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    .       0       0       0       0       0       0

Summary by RUN
==============

`diff` is the simulated `CP` from `nonmem` minus the simulated `CP` from `mrgsim`

``` r
group_by(comp,ID) %>% summarise(mean = mean(diff), max = max(diff), min = min(diff))
```

    . # A tibble: 20 x 4
    .       ID  mean   max   min
    .    <int> <dbl> <dbl> <dbl>
    .  1     1     0     0     0
    .  2     2     0     0     0
    .  3     3     0     0     0
    .  4     4     0     0     0
    .  5     5     0     0     0
    .  6     6     0     0     0
    .  7     7     0     0     0
    .  8     8     0     0     0
    .  9     9     0     0     0
    . 10    10     0     0     0
    . 11    11     0     0     0
    . 12    12     0     0     0
    . 13    13     0     0     0
    . 14    14     0     0     0
    . 15    15     0     0     0
    . 16    16     0     0     0
    . 17    17     0     0     0
    . 18    18     0     0     0
    . 19    19     0     0     0
    . 20    20     0     0     0

``` r
comp_plot <- function(comp) {
  ggplot(data = comp) + 
    ggtitle(paste0("ID: ", comp$ID[1]), subtitle="Line: mrgsolve, Point: NONMEM") + 
    geom_point(aes(time,NONMEM),color = "firebrick") + 
    geom_line(aes(time,MRGSIM,group = ID)) +
    theme_bw() + ylab("Simulated value") + 
    scale_x_continuous(breaks = seq(0,130,24))  
}


runs <- mutate(runs, plot = map(comp, comp_plot))
```

Results
=======

``` r
get_title <- function(i) unlist(slice(runs,i) %>% select(descr))
```

Bolus with additional
---------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl ID
    . 1    0   1 100    1 24    3  1
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-39-1.png)

Bolus with lag time and bioav
-----------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl  LAGT BIOAV ID
    . 1    0   2 100    1 24    3 12.13  2.23  2
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-40-1.png)

Infusion with additional
------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate ID
    . 1    0   2 100    1 24    3   10  3
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-41-1.png)

Infusion doses to depot, with additional
----------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl     rate ID
    . 1    0   1 100    1 24    3 8.333333  4
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-42-1.png)

Infusion doses, with additional and lag time
--------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl     rate ID
    . 1    0   1 100    1 24    3 8.333333  5
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-43-1.png)

Infusion doses, with lag time and bioav factor
----------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl     rate ID
    . 1    0   1 100    1 24    3 8.333333  6
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-44-1.png)

Infusion doses, with lag time and bioav factor
----------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate LAGT BIOAV ss ID
    . 1    0   2 100    1 24    3   10 3.16 0.412  1  7
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-45-1.png)

Infusion doses at steady-state, with lag time and bioav factor
--------------------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate BIOAV ss ID
    . 1    0   2 100    1 12    4    2 0.812  1  8
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-46-1.png)

Infusion doses, with lag time and bioav factor
----------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate ss ID
    . 1    0   2 100    1 12    3    2  1  9
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-47-1.png)

Infusion doses at steady state, II &lt; DUR, no bioav factor
------------------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl   rate ss ID
    . 1    0   2 100    1  6   12 8.3333  1 10
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-48-1.png)

Infusion doses at steady state where II == DUR, with bioav factor
-----------------------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate BIOAV ss ID
    . 1    0   2 100    1 10    8 4.12 0.412  1 11
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-49-1.png)

Infusion doses at steady state, where II == DUR
-----------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate ss ID
    . 1    0   2 100    1 10    8   10  1 12
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-50-1.png)

Bolus doses at steady state, with bioav factor and lag time
-----------------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl LAGT BIOAV ss ID
    . 1    0   2 100    1 24    3    4 0.412  1 13
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-51-1.png)

Bolus doses with lag time and bioavability factor
-------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl LAGT BIOAV ID
    . 1    0   2 100    1 24    3    5 0.412 14
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-52-1.png)

Bolus then infusion
-------------------

    . $ev
    . Events:
    .   time cmt amt evid LAGT ii addl rate ID
    . 1    0   2 100    1    1  0    0    0 15
    . 2   13   1  50    1    0 24    2   24 15
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-53-1.png)

Infusion with modeled duration, lag time, and bioav factor
----------------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid rate DUR2 MODE ii addl LAGT BIOAV ID
    . 1    0   2 100    1   -2    9    2 24    3    5  0.61 16
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-54-1.png)

Infusion with modeled duration, at steady state with bioav factor
-----------------------------------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid rate DUR2 MODE ii addl ss BIOAV ID
    . 1    0   2 100    1   -2    9    2 24    3  1  0.61 17
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-55-1.png)

Reset and dose (EVID 4) with additional
---------------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate BIOAV ID
    . 1    0   1 100    1 12    2   50  0.61 18
    . 2   50   1 120    4 12    3    0  0.50 18
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-56-1.png)

Reset (EVID 3) with additional
------------------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl rate BIOAV ID
    . 1    0   1 100    1 12    3   50  0.61 19
    . 2   50   2   0    3  0    0    0  0.00 19
    . 3   54   1 120    1 16    2    0  0.00 19
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-57-1.png)

Steady state 1 and 2
--------------------

    . $ev
    . Events:
    .   time cmt amt evid ii addl ss ID
    . 1    0   1 100    1 24    3  1 20
    . 2   12   1  50    1 24    3  2 20
    . 
    . $plot

![](img/nmtest4-unnamed-chunk-58-1.png)

Control stream
==============

``` r
writeLines(readLines("model/1001.ctl"))
```

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
       
       $TABLE FILE=TAB ID TIME EVID CP IPRED PRED DV NOPRINT ONEHEADER NOAPPEND
       
       $SIMULATION (2674474) ONLYSIMULATION

Session Info
============

``` r
options(width = 120)
devtools::session_info()
```

    . ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
    .  setting  value                       
    .  version  R version 3.3.3 (2017-03-06)
    .  os       Ubuntu 14.04.5 LTS          
    .  system   x86_64, linux-gnu           
    .  ui       X11                         
    .  language (EN)                        
    .  collate  en_US.UTF-8                 
    .  ctype    en_US.UTF-8                 
    .  tz       Etc/UTC                     
    .  date     2018-12-20                  
    . 
    . ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
    .  package       * version     date       lib source                                       
    .  assertthat      0.2.0       2017-04-11 [1] CRAN (R 3.3.3)                               
    .  backports       1.1.3       2018-12-14 [1] CRAN (R 3.3.3)                               
    .  bindr           0.1.1       2018-03-13 [1] CRAN (R 3.3.3)                               
    .  bindrcpp      * 0.2.2       2018-03-29 [1] CRAN (R 3.3.3)                               
    .  callr           3.1.0       2018-12-10 [1] CRAN (R 3.3.3)                               
    .  cli             1.0.1       2018-09-25 [1] CRAN (R 3.3.3)                               
    .  colorspace      1.3-2       2016-12-14 [4] CRAN (R 3.3.2)                               
    .  crayon          1.3.4       2017-09-16 [1] CRAN (R 3.3.3)                               
    .  desc            1.2.0       2018-05-01 [1] CRAN (R 3.3.3)                               
    .  devtools        2.0.1       2018-10-26 [1] CRAN (R 3.3.3)                               
    .  digest          0.6.18      2018-10-10 [1] CRAN (R 3.3.3)                               
    .  dplyr         * 0.7.8       2018-11-10 [1] CRAN (R 3.3.3)                               
    .  evaluate        0.12        2018-10-09 [1] CRAN (R 3.3.3)                               
    .  fansi           0.4.0       2018-10-05 [1] CRAN (R 3.3.3)                               
    .  fs              1.2.6       2018-08-23 [1] CRAN (R 3.3.3)                               
    .  ggplot2       * 3.1.0       2018-10-25 [1] CRAN (R 3.3.3)                               
    .  glue            1.3.0       2018-07-17 [1] CRAN (R 3.3.3)                               
    .  gtable          0.2.0       2016-02-26 [4] CRAN (R 3.3.2)                               
    .  highr           0.6         2016-05-09 [2] CRAN (R 3.3.2)                               
    .  hms             0.4.2       2018-03-10 [1] CRAN (R 3.3.3)                               
    .  htmltools       0.3.6       2017-04-28 [1] CRAN (R 3.3.3)                               
    .  knitr           1.21        2018-12-10 [1] CRAN (R 3.3.3)                               
    .  labeling        0.3         2014-08-23 [4] CRAN (R 3.3.2)                               
    .  lattice         0.20-34     2016-09-06 [4] CRAN (R 3.3.2)                               
    .  lazyeval        0.2.1       2017-10-29 [1] CRAN (R 3.3.3)                               
    .  magrittr        1.5         2014-11-22 [2] CRAN (R 3.3.2)                               
    .  MASS            7.3-51.1    2018-11-01 [1] CRAN (R 3.3.3)                               
    .  memoise         1.0.0       2016-01-29 [4] CRAN (R 3.3.2)                               
    .  metrumrg        5.57        2017-10-14 [1] Github (metrumresearchgroup/metrumrg@2e5a541)
    .  mrgsolve      * 0.8.12.9000 2018-12-19 [1] Github (metrumresearchgroup/mrgsolve@4286b09)
    .  munsell         0.5.0       2018-06-12 [1] CRAN (R 3.3.3)                               
    .  pillar          1.3.1       2018-12-15 [1] CRAN (R 3.3.3)                               
    .  pkgbuild        1.0.2       2018-10-16 [1] CRAN (R 3.3.3)                               
    .  pkgconfig       2.0.2       2018-08-16 [1] CRAN (R 3.3.3)                               
    .  pkgload         1.0.2       2018-10-29 [1] CRAN (R 3.3.3)                               
    .  plyr            1.8.4       2016-06-08 [4] CRAN (R 3.3.2)                               
    .  prettyunits     1.0.2       2015-07-13 [1] CRAN (R 3.3.3)                               
    .  processx        3.2.1       2018-12-05 [1] CRAN (R 3.3.3)                               
    .  ps              1.2.1       2018-11-06 [1] CRAN (R 3.3.3)                               
    .  purrr         * 0.2.5       2018-05-29 [1] CRAN (R 3.3.3)                               
    .  R6              2.3.0       2018-10-04 [1] CRAN (R 3.3.3)                               
    .  Rcpp            1.0.0       2018-11-07 [1] CRAN (R 3.3.3)                               
    .  RcppArmadillo   0.9.200.5.0 2018-11-28 [1] CRAN (R 3.3.3)                               
    .  readr         * 1.3.0       2018-12-11 [1] CRAN (R 3.3.3)                               
    .  remotes         2.0.2       2018-10-30 [1] CRAN (R 3.3.3)                               
    .  reshape         0.8.8       2018-10-23 [1] CRAN (R 3.3.3)                               
    .  rlang           0.3.0.1     2018-10-25 [1] CRAN (R 3.3.3)                               
    .  rmarkdown       1.11        2018-12-08 [1] CRAN (R 3.3.3)                               
    .  rprojroot       1.3-2       2018-01-03 [1] CRAN (R 3.3.3)                               
    .  rstudioapi      0.8         2018-10-02 [1] CRAN (R 3.3.3)                               
    .  scales          1.0.0       2018-08-09 [1] CRAN (R 3.3.3)                               
    .  sessioninfo     1.1.1       2018-11-05 [1] CRAN (R 3.3.3)                               
    .  stringi         1.2.4       2018-07-20 [1] CRAN (R 3.3.3)                               
    .  stringr         1.3.1       2018-05-10 [1] CRAN (R 3.3.3)                               
    .  sys             2.1         2018-11-13 [1] CRAN (R 3.3.3)                               
    .  testthat        2.0.1       2018-10-13 [1] CRAN (R 3.3.3)                               
    .  tibble          1.4.2       2018-01-22 [1] CRAN (R 3.3.3)                               
    .  tidyr         * 0.8.2       2018-10-28 [1] CRAN (R 3.3.3)                               
    .  tidyselect      0.2.5       2018-10-11 [1] CRAN (R 3.3.3)                               
    .  usethis         1.4.0       2018-08-14 [1] CRAN (R 3.3.3)                               
    .  utf8            1.1.4       2018-05-24 [1] CRAN (R 3.3.3)                               
    .  withr           2.1.2       2018-03-15 [1] CRAN (R 3.3.3)                               
    .  xfun            0.4         2018-10-23 [1] CRAN (R 3.3.3)                               
    .  XML             3.98-1.16   2018-08-19 [1] CRAN (R 3.3.3)                               
    .  yaml            2.2.0       2018-07-25 [1] CRAN (R 3.3.3)                               
    . 
    . [1] /data/Rlibs
    . [2] /usr/local/lib/R/site-library
    . [3] /usr/lib/R/site-library
    . [4] /usr/lib/R/library
