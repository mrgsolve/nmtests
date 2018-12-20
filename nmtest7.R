##' ---
##' title: "Tests with NONMEM"
##' author: "Metrum Research Group, LLC"
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
knitr::opts_chunk$set(comment = '.', fig.path = "img/nmtest4-")
#+

carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")

##' # Functions
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

nmtest_plot <- function(x) {
  plot(x, CP~time)  
}

##' 
##' 
##' # The `mrgsim` model

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
mod <- update(mod, end=130, delta = 1)

##'
##'
##'
##' # Assemble the scenarios
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
push_back(env, ev, "Bolus with lag time and bioav")
#+
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, cmt = 2) 
push_back(env,ev, "Infusion with additional")
#+
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/12, cmt = 1) 
push_back(env,ev,"Infusion doses to depot, with additional")
#+
ev4 <- ev(amt = 100, ii = 24, addl=3, rate = 100/10, LAGT = 4.15, cmt = 2) 
push_back(env,ev,"Infusion doses, with additional and lag time")
#+
ev5 <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.25, BIOAV = 0.412, cmt = 2) 
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
ev <- ev(amt = 100, ii = 24, addl=3,  LAGT = 4, BIOAV = 0.412, ss = 1, cmt = 2) 
push_back(env,ev,"Bolus doses at steady state, with bioav factor and lag time")
#+
ev <- ev(amt = 100, ii = 24, addl=3,  LAGT = 5, BIOAV = 0.412, cmt = 2) 
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
  ev(amt = 0, evid = 3, time = 50, cmt = 2) + 
  ev(amt = 120, ii = 16, addl = 2, time = 54)
push_back(env,ev,"Reset (EVID 3) with additional")
#+
ev <- 
  ev(amt = 100, ii = 24, addl = 3, ss = 1)  + 
  ev(amt = 50,  ii = 24, addl = 3, ss = 2, time = 12)
push_back(env,ev,"Steady state 1 and 2")
#+
update_id <- function(ev,id) mutate(ev, ID = id)

runs <- data_frame(ev = env$ev, descr = env$descr)
runs <- mutate(runs, ID = seq(n()))
runs <- mutate(runs,ev = map2(ev,ID, update_id))
runs <- mutate(runs, sims = mclapply(ev, sim, x = mod))

runs <- mutate(runs, data = map(sims, to_data_set))

#+
data <- runs[["data"]] %>% bind_rows()

sv(data, "data/1001.csv")

##' # Simulate with `nonmem`
out <- run(1001)

##' # Overall Summary
##' 
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


##' This is the `nonmem` minus `mrgsim` summary
summary(comp$diff)

##' # Summary by RUN
##' 
##' `diff` is the simulated `CP` from `nonmem` minus the simulated
##' `CP` from `mrgsim`
group_by(comp,ID) %>% summarise(mean = mean(diff), max = max(diff), min = min(diff))

comp_plot <- function(comp) {
  ggplot(data = comp) + 
    ggtitle(paste0("ID: ", comp$ID[1]), subtitle="Line: mrgsolve, Point: NONMEM") + 
    geom_point(aes(time,NONMEM),color = "firebrick") + 
    geom_line(aes(time,MRGSIM,group = ID)) +
    theme_bw() + ylab("Simulated value") + 
    scale_x_continuous(breaks = seq(0,130,24))  
}


runs <- mutate(runs, plot = map(comp, comp_plot))

##' # Results

get_title <- function(i) unlist(slice(runs,i) %>% select(descr))

##' ## `r i <- 1; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 2; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 3; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 4; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 5; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 6; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 7; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 8; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 9; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 10; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 11; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 12; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 13; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 14; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 15; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 16; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 17; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 18; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 19; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()

##' ## `r i <- 20; get_title(i)`
#+ echo = FALSE
slice(runs, i) %>% select(ev,plot) %>% flatten()


##' # Control stream
#+ comment = "  "
writeLines(readLines("model/1001.ctl"))

##'
##' 
##' # Session Info
#+ 
options(width = 120)
devtools::session_info()

