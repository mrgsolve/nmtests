#' ---
#' title: "Benchmark test with mrgsolve and NONMEM"
#' author: "Metrum Research Group"
#' date: ""
#' output: 
#'   github_document:
#'     toc: true
#'     toc_depth: 2
#'   pdf_document:
#'     toc: true
#'     toc_depth: 2
#'     number_sections: true
#'     latex_engine: xelatex
#'     extra_dependencies:
#'       fontenc: T1
#'       mathdesign: utopia
#' ---

#' # Introduction
#' 
#' This document runs simulations from a pharmacokinetic model using both 
#' NONMEM and mrgsolve and compares the results. The benchmarks in this 
#' set focus on dosing events (bolus and infusion), bioavailability, 
#' lag times, reset and steady state, 
#' 
#' All of the relevant code is presented so that the user can trace how 
#' the simulations are performed.  The complete source code can be viewed
#' [here](dosing/dosing-vignette.R).
#' 
#' The bottom line results are presented in graphical
#' form  [here](#results) and numeric form [here](#numeric-summary).
#' 
#' 
#' # mrgsolve package version
packageVersion("mrgsolve")
#' 
#' # Setup

Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")

#+ message = FALSE
mrgsolve.loc <- NULL # "/data/home/Rlibs/"
library(mrgsolve, lib.loc = mrgsolve.loc)
# --------------------------------
#+ message = FALSE
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(tidyr)
library(jsonlite)
library(tools)
library(here)
library(knitr)

#+ echo = FALSE
knitr::opts_chunk$set(comment = '.', fig.path = "results/img/dosing-")
#+

stopifnot(grepl("PsN", system("execute --version", intern = TRUE)))
stopifnot(file.exists("locate-dosing"))
source(here("shared", "tools.R"))
source(here("shared", "data.R"))

#' # The `mrgsim()` model
mod <- mread_cache("1001.mod")
mod <- update(mod, end = 130, delta = 1)

runs <- data_1001()
runs <- mutate(runs, sims = lapply(ev, sim, x = mod))
runs <- mutate(runs, data = map(sims, to_data_set))

#+
data <- bind_rows(runs$data)
data$CP <- NULL

fsave(data, "data/1001.csv")

#' # Simulate with `nonmem`
nm <- psn_execute(1001)

#' # Numeric Summary
#' 
#' Look at the difference between simulated values from mrgsolve and NONMEM.
#' 
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

#' ## Overall
#' 
#' This is the `nonmem` minus `mrgsim()` summary
summary(comp$diff)

#' ## Summary by scenario number
#' 
#' `diff` is the simulated `CP` from `nonmem` minus the simulated `CP` 
#' from `mrgsim()`
group_by(comp, ID) %>% 
  summarise(mean = mean(diff), max = max(diff), min = min(diff)) %>% 
  kable()

runs <- mutate(runs, plot = map(comp, comp_plot))

#+ include = FALSE
sims <- pull(runs, sims) %>% lapply(as_tibble) %>% rbindlist()
fwrite(x = comp, file = "results/1001.csv")
fwrite(x = sims, file = "results/1001-sims.csv")
fwrite(x = nm, file = "results/1001-nm.csv")
run_key <- distinct(runs, ID, descr) %>% mutate(descr = unlist(descr))
run_key <- select(run_key, ID, descr)
fwrite(x = run_key, file = "results/1001-run-key.csv")

meta <- list()
meta$data <- list(
  file = "data/1001.csv", 
  md5 = md5sum("data/1001.csv")
)
meta$ctl <- list(
  file = "1001.ctl", 
  md5 = md5sum("1001.ctl")
)
meta$mod <- list(
  file  = "1001.mod", 
  md5 = md5sum("1001.mod")
)
meta$mrgsolve <- list(
  file = "results/1001-sims.csv", 
  md5 = md5sum("results/1001-sims.csv")
)
meta$nonmem <- list(
  file = "results/1001-nm.csv", 
  md5 = md5sum("results/1001-nm.csv")
)
meta$compare <- list(
  file = "results/1001.csv", 
  md5 = md5sum("results/1001.csv")
)
meta$key <- list(
  file = "results/1001-run-key.csv", 
  md5 = md5sum("results/1001-run-key.csv")
)
write_json(
  x = meta, 
  path = "results/1001.json",
  pretty = TRUE
)
#+


#' # Results
#+ echo = FALSE
get_title <- function(i) {
  de <- unlist(slice(runs,i) %>% select(descr))
  paste0(i, ": ", de)
}

#' ## `r i <- 1; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 2; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 2; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 4; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 5; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 7; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 8; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 9; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 10; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 11; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 12; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 13; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 14; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 15; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 16; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 17; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 18; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 19; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 20; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 21; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 22; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 23; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)

#' ## `r i <- 24; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% map(1)


#+ include = FALSE
if(i != nrow(runs)) stop("missing output")


#' # Control stream
#+ comment = "  "
writeLines(readLines("1001/1001.lst"))

#' # Session Info
#+ 
options(width = 120)
sessionInfo()
