##' ---
##' title: "Bemchmark test with mrgsolve and NONMEM"
##' author: "Metrum Research Group"
##' date: ""
##' output: 
##'   github_document:
##'     toc: true
##'     toc_depth: 2
##'   pdf_document:
##'     toc: true
##'     toc_depth: 2
##'     number_sections: true
##'     latex_engine: xelatex
##'     extra_dependencies:
##'       fontenc: T1
##'       mathdesign: utopia
##'     
##' ---

##' # Introduction
##' 
##' This document runs simulations from a pharmacokinetic model using both 
##' NONMEM and mrgsolve and compares the results.  
##' 
##' All of the relevant code is presented so that the user can trace how 
##' the simulations are performed.  The complete source code can be viewed
##' [here](nmtest7.R).
##' 
##' The bottom line results are presented in graphical
##' form  [here](#results) and numeric form [here](#numeric-summary).
##' 
##' 
##' 
##' # Setup

Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")

#+ message = FALSE
.libPaths("/data/Rlibs")
library(mrgsolve)
library(dplyr)
library(readr)
library(ggplot2)
library(parallel)
library(purrr)
library(tidyr)







carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")

##' # Functions
##' 
##' These functions assemble data sets, run simulations, and 
##' gather outputs.  All scenarios are handled in exactly the 
##' same way.
##' 
##' 
##' ## Save `mrgsim` output as a `nonmem` input data set
to_data_set <- function(x,id=NULL) {
  x <- as.data.frame(x)
  x <- mutate(x, C = '.', DV = '.', cmt = if_else(cmt==0, 2, cmt))
  x <- dplyr::select(x, "C", everything())
  if(is.numeric(id)) x <- mutate(x,ID = id)
  x
}

##' ## Save the `nonmem` input data set
sv <- function(x,file) {
  write.csv(file = file, row.names = FALSE, quote = FALSE, x)
}

##' ## Run `nonmem`
run <- function(number) {
  metrumrg::NONR(number, project = "model", 
                 command = "/opt/NONMEM/nm74/nmqual/autolog.pl", 
                 checkrunno=FALSE)
  return(tabread(number))
}

##' ## Read in `nonmem` simulation results
tabread <- function(number) {
  tab <- file.path("model", number, "TAB")
  if(file.exists(tab)) return(read_table(tab, skip=1))
  stop("the run failed")
}

##' ## Simulate a scenario with `mrsim`
sim <- function(x, e,...) {
  mrgsim(x, events = e, carry.out = carry, digits = 5, recsort=3, ...) 
}

##' 
##' # The `mrgsim` model

# code <- '
# $SET req = ""
# $PARAM CL = 1.1, V = 20, KA = 1.5
# LAGT = 0, MODE = 0, DUR2 = 2, RAT2 = 10, BIOAV = 1
# 
# $PKMODEL cmt = "GUT CENT", depot = TRUE
# 
# $MAIN
# 
# F_CENT = BIOAV;
# ALAG_CENT = LAGT;
# 
# if(MODE==1) R_CENT = RAT2;
# if(MODE==2) D_CENT = DUR2;
# 
# 
# $TABLE
# capture DV = (CENT/(V/1000));
# capture CP = DV;
# 
# $CAPTURE LAGT MODE DUR2 RAT2 BIOAV
# '
# 
# mod <- mcode_cache("tests1", code)
# mod <- update(mod, end=130, delta = 1)

##'
##'
##'
##' # Assemble the scenarios
##' 
##' There is a lot of code here.  See the [results](#results) section to 
##' see input data objects next to simulated data output from 
##' mrgsolve and NONMEM.
##' 
##' - Doses into `cmt` 2 are intravascular and doses into `cmt` 1 are extravascular
##' - `LAGT` sets the dosing lag time
##' - `BIOAV` sets the bioavailability fraction
##' 
##'
##'
env <- new.env()
env$ev <- list()
env$descr <- list()
push_back <- function(env,ev,descr) {
  n <- length(env$ev)+1
  m <- length(env$descr)+1
  env$ev[[n]] <- ev
  env$descr[[m]] <- descr
}

#+ 
ev <- ev(amt = 100, ii = 24, addl = 3) 
push_back(env,ev, "Bolus with additional")

#+
ev <- ev(amt = 100, ii = 24, addl = 3, LAGT = 12.13, BIOAV = 2.23, cmt = 2) 
push_back(env, ev,"Bolus with lag time and bioav")
#+
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, cmt = 2) 
push_back(env,ev,"Infusion with additional")
#+
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/12, cmt = 1) 
push_back(env,ev,"Infusion doses to depot, with additional")
#+
ev <- ev(amt = 100, ii = 24, addl=3, rate = 100/10, LAGT = 4.15, cmt = 2) 
push_back(env,ev,"Infusion doses, with additional and lag time")
#+
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.25, BIOAV = 0.412, cmt = 2) 
push_back(env,ev,"Infusion doses, with lag time and bioav factor")
#+
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.16, BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses, with lag time and bioav factor")
#+
ev <- ev(amt = 100, ii = 12, addl = 4, rate = 100/50, BIOAV = 0.812, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady-state, with lag time and bioav factor")
#+
ev <- ev(amt = 100, ii = 12, addl = 3, rate = 100/50, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses, with lag time and bioav factor")
#+
ev <- ev(amt = 100, ii = 6, addl = 12, rate = signif(100/12,5), ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady state, II < DUR, no bioav factor")
#+
ev <- ev(amt = 100, ii = 10, addl = 8, rate = 0.412*100/10,  BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady state where II == DUR, with bioav factor")
#+
ev <- ev(amt = 100, ii = 10, addl = 8, rate = 100/10, ss = 1, cmt = 2) 
push_back(env,ev,"Infusion doses at steady state, where II == DUR")
#+
ev <- ev(amt = 100, ii = 24, addl = 3,  LAGT = 4, BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Bolus doses at steady state, with bioav factor and lag time")
#+
ev <- ev(amt = 100, ii = 24, addl = 3,  LAGT = 5, BIOAV = 0.412, cmt = 2) 
push_back(env,ev,"Bolus doses with lag time and bioavability factor")
#+
ev <- 
  ev(amt = 100, cmt = 2, LAGT = 1) + 
  ev(time = 13, amt = 50, ii = 24, addl = 2, rate = 24)
push_back(env,ev,"Bolus then infusion")
#+
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, LAGT = 5, BIOAV = 0.61)
push_back(env,ev,"Infusion with modeled duration, lag time, and bioav factor")
#+
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, ss = 1, BIOAV = 0.61)
push_back(env,ev,"Infusion with modeled duration, at steady state with bioav factor")
#+
ev <- 
  ev(amt = 100, ii = 12, addl = 2, rate = 50, BIOAV = 0.61) + 
  ev(amt = 120, evid = 4, time = 50, BIOAV = 0.5, ii = 12, addl = 3)
push_back(env,ev,"Reset and dose (EVID 4) with additional")
#+
ev <- 
  ev(amt = 100, ii = 12, addl = 3, rate = 50, BIOAV = 0.61) + 
  ev(amt = 0, evid = 3, time = 50, cmt = 2, BIOAV=1) + 
  ev(amt = 120, ii = 16, addl = 2, time = 54, BIOAV=1)
push_back(env,ev,"Reset (EVID 3) with additional")

#+
ev <- 
  ev(amt = 100, ii = 24, addl = 3, ss = 1)  + 
  ev(amt = 50,  ii = 24, addl = 3, ss = 2, time = 12)
push_back(env,ev,"Steady state 1 and 2")
#+ 
ev <- ev(amt = 0, rate = 100,  ss=1)
push_back(env,ev,"Steady state infusion")

#+
update_id <- function(ev,id) mutate(as.data.frame(ev), ID = id)
runs <- tibble(ev = env$ev, descr = env$descr)
runs <- mutate(runs, ID = seq(n()))
runs <- mutate(runs,ev = map2(ev,ID, update_id))

evs <- runs[["ev"]]
lapply(evs[15], function(x) {
  print(x)
  expand_observations(x,obs)
})

dd <- evs[15]
expand_observations(dd[[1]],obs)
mod <- mrgsolve:::house()
mrgsim(mod,dd[[1]])


x <- ev(time=c(0,13), amt = c(100,50), cmt=c(2,1), evid=1, LAGT = c(1,0), rate = c(0,24), addl = c(0,2), ii=c(0,24), ID=15)
y <- evs[[15]]

x
expand_observations(as.data.frame(x),obs)
mrgsim_d(mod,x)

obs <- seq(1,130)
runs <- mutate(runs, ev = map(ev, expand_observations, times = obs))



runs$ev[[15]]
x <- runs$ev[[16]]
expand_observations(x,obs)


dat <- expand_observations(x,obs)




runs <- mutate(runs, sims = mclapply(ev, sim, x = mod))

runs <- mutate(runs, data = map(sims, to_data_set))

#+
data <- runs[["data"]] %>% bind_rows()


sv(data, "data/1001.csv")

##' # Simulate with `nonmem`
out <- run(1001)
