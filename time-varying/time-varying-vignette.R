#' ---
#' title: "Time-varying covariates with mrgsolve and NONMEM"
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
#'     
#' ---

#' # Introduction
#' 
#' This document runs simulations from a pharmacokinetic model that 
#' involves time-varying covariates and compares the result.  A 
#' more-comprehensive comparision of different dosing scenarios is 
#' provided in this repository, but in another document.
#' 
#' All of the relevant code is presented so that the user can trace how 
#' the simulations are performed.  The complete source code can be viewed
#' [here](nmtest_time_varying.R).
#' 
#' 
#' 
#' # Setup

Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")

#+ message = FALSE
library(mrgsolve)
library(dplyr)
library(data.table)
library(ggplot2)
library(parallel)
library(purrr)
library(here)
library(tools)
library(jsonlite)

source(here("shared", "tools.R"))

#+ echo = FALSE
knitr::opts_chunk$set(comment = '.', fig.path = "results/img/time-varying-")
#+

#' 
#' # The `mrgsim` model
#' 
#' - `WT` is included as a covariate on `CL`
#' - Note that we capture both `WT` and `CL` for comparison later
mod <- mread("2001.mod")
mod <- update(mod, end = 130, delta = 4)

#' # Assemble the scenarios
#' 
#' - We are only testing one dosing intervention here.  The focus 
#' is seeing what happens with time-varying covariates. 
#' 
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

update_id <- function(ev,id) mutate(ev, ID = id)

runs <- tibble(ev = env$ev, descr = env$descr)
runs <- mutate(runs, ID = seq(n()))
runs <- mutate(runs,ev = map2(ev,ID, update_id))
runs <- mutate(runs, sims = lapply(ev, sim, x = mod))
runs <- mutate(runs, data = map(sims, to_data_set))


#' # The input data set
#' 
#+
data <- bind_rows(runs$data)

#+
head(data)

#' Add a weight column to the data
set.seed(10010)
wt <- distinct(data, ID, time) %>% 
  mutate(WT = round(runif(n(), 60,80),2))

data <- select(data, -WT) %>% left_join(wt)
data <- select(data, C, ID, time, WT, everything())

#+
head(data)

#' 
#' And resimulate so that `CL` is a function of `WT`
#' 
sims <- mrgsim_d(mod, data, digits = 5, recsort = 4, output = "df")

#' Notice that we simulated with `recsort` equal to 4.  This is important
#' to control record sort order when scheduling doses via `addl`.  This is 
#' important even when simulating from a data set that includes observation 
#' records. 
#' 

fsave(data, "data/2001.csv")

#' # Simulate with `nonmem`
nm <- psn_execute(2001)

#' # Summary
comp <- mutate(
  sims, 
  NONMEM = nm$CP, 
  MRGSIM = CP, 
  CLNM = nm$CL, 
  CLMG = CL, 
  WTNM = nm$WT, 
  WTMG = WT
)

#+ 
head(comp)

#' ## Numerical summary
summary(comp$NONMEM - comp$MRGSIM)

#' ## Concentration versus time
#+ 
ggplot(comp) + 
  geom_point(aes(time, NONMEM), col = "firebrick", size = 3) + 
  geom_line(aes(time, MRGSIM), lwd=1) + theme_bw() + 
  ylab("Concentration")

#' ## Clearance versus time
#+  
ggplot(comp) + 
  geom_point(aes(time, CLNM), col = "firebrick", size = 3) + 
  geom_line(aes(time, CLMG), lwd=1) + theme_bw() + 
  ylab("Clearance") + 
  ylim(0, NA)

#' ## Weight versus time
#+ 
ggplot(comp) + 
  geom_point(aes(time, WTNM),col = "firebrick", size = 3) + 
  geom_line(aes(time, WTMG), lwd = 1) + theme_bw() + 
  ylim(0, NA) + 
  ylab("Weight")

#+ include = FALSE
fwrite(x = comp, file = "results/2001.csv")
fwrite(x = sims, file = "results/2001-sims.csv")
fwrite(x = nm, file = "results/2001-nm.csv")

meta <- list()
meta$md5 <- list()
meta$data <- list(
  file = "data/2001.csv", 
  md5 = md5sum("data/2001.csv")
)
meta$ctl <- list(
  file = "2001.ctl", 
  md5 = md5sum("2001.ctl")
)
meta$mod <- list(
  file  = "2001.mod", 
  md5 = md5sum("2001.mod")
)
meta$mrgsolve <- list(
  file = "results/2001-sims.csv", 
  md5 = md5sum("results/2001-sims.csv")
)
meta$nonmem <- list(
  file = "results/2001-nm.csv", 
  md5 = md5sum("results/2001-nm.csv")
)
meta$compare <- list(
  file = "results/2001.csv", 
  md5 = md5sum("results/2001.csv")
)
meta$date <- date()
meta$user <- Sys.info()[["user"]]
write_json(
  x = meta, 
  path = "results/2001.json",
  pretty = TRUE
)
#+

#' # Control stream
#+ comment = "  "
writeLines(readLines("2001/2001.lst"))

#' # Session Info
#+ 
options(width = 120)
sessionInfo()
