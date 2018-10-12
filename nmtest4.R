##' ---
##' title: "Tests with NONMEM"
##' author: "Metrum Research Group, LLC"
##' date: ""
##' output: 
##'   github_document:
##'     toc: true
##'   pdf_document:
##'     toc: true
##'     number_sections: true
##'     latex_engine: xelatex
##'     extra_dependencies:
##'       fontenc: T1
##'       mathdesign: utopia

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
ev <- ev(amt = 100, ii = 24, addl = 3) 
ev
out1 <- sim(mod,ev)
plot(out1)
data1 <- to_data_set(out1, 1)

##' ### Bolus doses, lag time and bioav factor
ev <- ev(amt = 100, ii = 24, addl = 3, LAGT = 12.13, BIOAV = 2.23, cmt = 2) 
ev
out1.1 <- sim(mod,ev)
plot(out1.1)
data1.1 <- to_data_set(out1.1, 1.1)

##' ### Infusion doses, with additional
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, cmt = 2) 
ev
out2 <- sim(mod,ev)
plot(out2)
data2 <- to_data_set(out2, 2)

##' ### Infusion doses to depot, with additional
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/12, cmt = 1) 
ev
out2.1 <- sim(mod,ev)
plot(out2.1)
data2.1 <- to_data_set(out2.1, 2.1)


##' ### Infusion doses, with additional and lag time
ev <- ev(amt = 100, ii = 24, addl=3, rate = 100/10, LAGT = 4.15, cmt = 2) 
ev
out3 <- sim(mod,ev)
plot(out3)
data3 <- to_data_set(out3, 3)

##' ### Infusion doses, with lag time and bioav factor
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.25, BIOAV = 0.412, cmt = 2) 
ev
out4 <- sim(mod,ev)
plot(out4)
data4 <- to_data_set(out4, 4)

##' ### Infusion doses at steady-state, with lag time and bioav factor
ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.16, BIOAV = 0.412, ss = 1, cmt = 2) 
ev
out5 <- sim(mod,ev)
plot(out5)
data5 <- to_data_set(out5, 5)

##' ### Infusion doses at steady state, II < DUR, with bioav factor
ev <- ev(amt = 100, ii = 6, addl = 12, rate = 100/10, BIOAV = 0.812, ss = 1, cmt = 2) 
ev
out6 <- sim(mod,ev)
plot(out6)
data6 <- to_data_set(out6, 6)

##' ### Infusion doses at steady state, II < DUR, no bioav factor
ev <- ev(amt = 100, ii = 6, addl = 12, rate = 100/10, ss = 1, cmt = 2) 
ev
out6.1 <- sim(mod,ev)
plot(out6.1)
data6.1 <- to_data_set(out6.1, 6.1)

##' ### Infusion doses at steady state where II is a multiple of DUR
ev <- ev(amt = 100, ii = 6, addl = 12, rate = signif(100/12,5), ss = 1, cmt = 2) 
ev
out6.2 <- sim(mod,ev)
plot(out6.2)
data6.2 <- to_data_set(out6.2, 6.2)

##' ### Infusion doses at steady state where II == DUR, with bioav factor
ev <- ev(amt = 100, ii = 10, addl = 8, rate = 100/10, LAGT = 0, BIOAV = 0.412, ss = 1, cmt = 2) 
ev
out7 <- sim(mod,ev)
plot(out7)
data7 <- to_data_set(out7, 7)

##' ### Infusion doses at steady state, where II == DUR
ev <- ev(amt = 100, ii = 10, addl = 8, rate = 100/10, ss = 1, cmt = 2) 
ev
out7.1 <- sim(mod,ev)
plot(out7.1)
data7.1 <- to_data_set(out7.1, 7.1)

##' ### Bolus doses at steady state, with bioav factor and lag time
ev <- ev(amt = 100, ii = 24, addl=3,  LAGT = 4, BIOAV = 0.412, ss = 1, cmt = 2) 
ev
out8 <- sim(mod,ev)
plot(out8)
data8 <- to_data_set(out8, 8)

##' ### Bolus doses with lag time and bioavability factor
ev <- ev(amt = 100, ii = 24, addl=3,  LAGT = 5, BIOAV = 0.412, cmt = 2) 
ev
out9 <- sim(mod,ev)
plot(out9)
data9 <- to_data_set(out9, 9)

##' ### Bolus / infusion
ev <- ev(amt = 100, cmt = 2, LAGT = 1) + ev(time = 13, amt = 50, ii = 24, addl = 2, rate = 24)
ev
out10 <- sim(mod,ev)
plot(out10)
data10 <- to_data_set(out10, 10)

##' ### Infusion with modeled duration, lag time, and bioav factor
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, LAGT = 5, BIOAV = 0.61)
ev
out11 <- sim(mod,ev)
plot(out11)
data11 <- to_data_set(out11,11)

##' ### Infusion with modeled duration, at steady state with bioav factor
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 12, addl = 5, ss = 1, BIOAV = 0.61)
ev
out12 <- sim(mod,ev)
plot(out12)
data12 <- to_data_set(out12,12)

##' ### Reset and dose (EVID 4) with additional
##' 
##' 
ev <- 
  ev(amt = 100, ii = 12, addl = 5, rate = 50, BIOAV = 0.61) + 
  ev(amt = 120, evid = 4, time = 80, BIOAV = 0.5, ii = 12, addl = 2)
ev
out13 <- sim(mod,ev)
plot(out13)
data13 <- to_data_set(out13,13)

##' ### Reset (EVID 3) with additional
##' 
##' 
ev <- 
  ev(amt = 100, ii = 12, addl = 3, rate = 50, BIOAV = 0.61) + 
  ev(amt = 0, evid = 3, time = 50, cmt = 2) + 
  ev(amt = 120, ii = 24, addl = 2, time = 54)
ev
out14 <- sim(mod,ev)
plot(out14)
data14 <- to_data_set(out14,14)

##' ### Steady state 1 and 2
##' 
##' 
ev <- 
  ev(amt = 100, ii = 24, addl = 3, ss = 1)  + 
  ev(amt = 50,  ii = 24, addl = 3, ss = 2, time = 12)
ev

out15 <- sim(mod,ev)
plot(out15)
data15 <- to_data_set(out15,15)

##' # Collect `mrgsim` output
sims <- list(out1,out1.1,out2,out2.1,out3,out4,out5,out6,out6.1,out6.2,out7,out7.1,
             out8,out9,out10,out11,out12,out13,out14,out15)
sims <- lapply(sims, as.data.frame)
sims <- bind_rows(sims)

##' # Create a single data set for `nonmem`
data <- bind_rows(data1,data1.1,data2,data2.1,data3,data4,data5,data6,data6.1,data6.2,data7,data7.1,
                  data8,data9,data10,data11,data12,data13,data14,data15)

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


##' # Summary by RUN
##' 
##' `diff` is the simulated `CP` from `nonmem` minus the simulated
##' `CP` from `mrgsim`
group_by(data,ID) %>% 
  mutate(diff = NM - MRGSIM) %>%
  summarise(mean = mean(diff), max = max(diff), min = min(diff))

##' # Plot
#+ fig.height = 10
ggplot(data = data) + 
  geom_point(aes(time,NM),color = "firebrick") + 
  geom_line(aes(time,MRGSIM,group = ID)) +
  facet_wrap(~ID, scales = "free_y", ncol = 2,strip.position = "right") + 
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



