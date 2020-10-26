##' ---
##' title: "Tests with NONMEM"
##' author: "Metrum Research Group, LLC"
##' date: ""
##' output: 
##'   pdf_document:
##'     toc: true
##'     number_sections: true
##'     latex_engine: xelatex
##'     extra_dependencies:
##'       fontenc: T1
##'       mathdesign: utopia
##'   github_document:
##'     toc: true
##'     
##' ---

##' \newpage

Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")

#+ message = FALSE
.libPaths("/data/Rlibs")
library(mrgsolve)
library(dplyr)
library(readr)
library(ggplot2)

#+ echo = FALSE
knitr::opts_chunk$set(comment = '.', fig.path = "img/nmtest4-")
#+

carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")

##' # Functions
##' 
##' 
##' ## Save `mrgsim` output as a `nonmem` input data set
to_data_set <- function(x,id = NULL) {
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
                 command = "/opt/NONMEM/nm73/nmqual/autolog.pl", 
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

##' \newpage
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
##'\newpage
##'
##' # Scenarios
##' 
##' ### Bolus doses, with additional
ev <- ev(amt = 100, ii = 24, addl = 3, ss = 1)  + ev(amt = 40, ii = 24, addl = 3, ss = 2, time = 12)
ev

out1 <- sim(mod,ev)
plot(out1)
data1 <- to_data_set(out1, 1)

##' # Collect `mrgsim` output
sims <- list(out1)
sims <- lapply(sims, as.data.frame)
sims <- bind_rows(sims)

##' # Create a single data set for `nonmem`
data <- bind_rows(data1)

sv(data, "data/1001.csv")

##' # Simulate with `nonmem`
out <- run(1001)


##' # Overall Summary
##' 
##' Dimensions for mrgsim and nonmem output
dim(out)
dim(sims)

##' This is the `nonmem` minus `mrgsim` summary
summary(out$CP - sims$CP)

data$NM <- out$CP
data$MRGSIM <- sims$CP

##' 
##' \newpage
##' 
##' # Summary by RUN
##' 
##' `diff` is the simulated `CP` from `nonmem` minus the simulated
##' `CP` from `mrgsim`
group_by(data,ID) %>% 
  mutate(diff = NM - MRGSIM) %>%
  summarise(mean = mean(diff), max = max(diff), min = min(diff))

##' # Plot
#+ fig.height = 8
ggplot(data = data) + 
  geom_point(aes(time,NM),color = "firebrick") + 
  geom_line(aes(time,MRGSIM,group = ID)) +
  facet_wrap(~ID, scales = "free_y", ncol = 3) + 
  theme_bw()



##' \newpage
##' 
##' # Control stream
#+ comment = "  "
writeLines(readLines("model/1001.ctl"))

##' \newpage
##' 
##' # Session Info
#+
devtools::session_info()



