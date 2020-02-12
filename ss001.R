

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

#+ echo = FALSE
knitr::opts_chunk$set(comment = '.', fig.path = "img/atol-")
#+

carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss", "WT")

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
  mrgsim(x, events = e, carry.out = carry, digits = 5, ...) 
}

##' 
##' # The `mrgsim` model
##' 
##' - `WT` is included as a covariate on `CL`
##' - Note that we capture both `WT` and `CL` for comparison later
##' 
##' 

code <- '
$SET req = ""
$PARAM TVCL = 1.1, V = 20, KA = 1.5
LAGT = 0, MODE = 0, DUR2 = 2, RAT2 = 10, BIOAV = 1, 
WT = 70

$PKMODEL cmt = "GUT CENT", depot = TRUE

$MAIN
double CL = TVCL*pow(WT/70,0.75); 

$TABLE
capture DV = (CENT/(V/1000));
capture CP = DV;

$CAPTURE LAGT MODE DUR2 RAT2 BIOAV WT CL
'

mod <- mcode_cache("tests_atol", code)
mod <- update(mod, end=148, delta = 0.5)

##'
##'
##'
##' # Assemble the scenarios
##' 
##' - We are only testing one dosing intervention here.  The focus 
##' is seeing what happens with time-varying covariates. 
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
ev <- ev(amt = 100, ii = 12, addl = 3) 
push_back(env,ev, "Bolus with additional")

update_id <- function(ev,id) mutate(ev, ID = id)

runs <- tibble(ev = env$ev, descr = env$descr)
runs <- mutate(runs, ID = seq(n()))
runs <- mutate(runs,ev = map2(ev,ID, update_id))
runs <- mutate(runs, sims = mclapply(ev, sim, x = mod))

runs <- mutate(runs, data = map(sims, to_data_set))


#' # The input data set
#' 
#' 
#+
data <- runs[["data"]] %>% bind_rows()

#+
head(data)


#' Add a weight column to the data
set.seed(10010)
wt <- distinct(data,ID,time) %>% mutate(WT = round(runif(n(), 60,80),2))

data <- data %>% mutate(WT = 70)

data <- select(data, C,ID,time,WT,everything())

#' 
#' And resimulate so that `CL` is a function of `WT`
#' 
sims <- mrgsim_d(mod, data, digits=5, recsort = 4)

#' Notice that we simulated with `recsort` equal to 4.  This is important
#' to control record sort order when scheduling doses via `addl`.  This is 
#' important even when simulating from a data set that includes observation 
#' records. 
#' 

data <- mutate(data, rate  = ifelse(evid != 0, 100/14, 0), cmt = 2, LAGT = 3)
data <- mutate(data, ss = ifelse(evid !=0, 1, 0), DUR2 = 0, RAT2 = 0)


sv(data, "data/ss001.csv")

##' # Simulate with `nonmem`
out <- run("ss001")



#' ## Concentration versus time
#+ 
ggplot(out) + 
  geom_line(aes(TIME,CP),col="firebrick",lwd=1) + 
  ylab("Concentration") + geom_vline(xintercept = c(12,14,18)) + 
  scale_x_continuous(breaks = seq(0,72,4))


 