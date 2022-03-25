##' ---
##' title: "Predictive check with mrgsolve and NONMEM"
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

#+ echo = FALSE
Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")

#+ message = FALSE
.libPaths("/data/Rlibs")
library(mrgsolve)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(data.table)

#+ echo = FALSE
knitr::opts_chunk$set(comment = '.', fig.path = "img/nmtest5-")
#+
##' # Helper functions
##' 
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
  if(file.exists(tab)) return(fread(tab, skip=1))
  stop("the run failed")
}


##' # Generate a template data set
e1 <- ev(amt = 100)  %>% ev_rep(ID = 1:30)
e2 <- ev(amt = 300)  %>% ev_rep(ID = 1:30)
data <- as_data_set(e1,e2)

set.seed(11010)
data <- mutate(data, WT = runif(n(), 50,95), DOSE = amt)

head(data)
count(data,DOSE)

##' ## The `mrgsolve` model
##' 
##' We'll use this to both simulate the predcheck for mrgsolve as well as 
##' to generate the data template
##' 
code <- '
$SET req = ""

$PARAM WT = 70, DOSE = 1

$PKMODEL cmt = "GUT CENT", depot = TRUE

$THETA
1 
20 
1.5

$OMEGA
0.09 
0.2 
0.5

$SIGMA 
0.02

$MAIN
double TVCL = THETA1*pow(WT/70,0.75);
double CL = TVCL*exp(ETA(1));

double TVV = THETA2*(WT/70);
double V  = TVV*exp(ETA(2));

double TVKA = THETA3;
double KA = TVKA*exp(ETA(3));

$TABLE
capture DV = (CENT/(V/1000))*exp(EPS(1));
capture CP = DV;
'

mod <- mcode_cache("tests5", code)


##' ## Sampling times
tg <- tgrid(0,24,4,add = c(0.5,1.5))

out <- mrgsim_d(
  mod, data, tgrid = tg, 
  carry.out = "WT,DOSE,evid,amt,cmt,time", Req = "DV"
) %>% as_data_frame()

##' Just one sample data set

head(out)

ggplot(out, aes(time,DV,col = factor(DOSE),group = ID)) + 
  geom_point() + geom_line() + theme(legend.position = "top")

##' Write this out to file
temp <- 
  out %>% 
  mutate(C=NA_character_,DV=NA_character_) %>% 
  select(C,everything())

write.csv(file = "data/104.csv", temp, quote = FALSE, row.names=FALSE, na = '.')

##' # Simulate
##' 
##' ## With NONMEM
foo <- run(104) %>% rename(IREP = V1, ID = V2, TIME = V3, EVID = V4, CP = V5, IPRED = V6, 
                           PRED = V7, DV = V8, DOSE = V9)

head(foo)

##' ## With mrgsolve
##' 
##' 1000 replicates
out <- parallel::mclapply(1:1000, function(i) {
  mod %>% mrgsim_d(temp,carry.out = "DOSE") %>% mutate(IREP = i,TIME = time)
}) %>% bind_rows()

head(out)


##' # Summarize
##' 
##' - Dose-normalize the concentrations (`DVN`)
##' - Take the 5th, 50th, and 95th percentiles
##' - Then take the median of the 5/50/95th percentiles
##' 
##' 
##' ## The NONMEM data
##' 
sum_nm <- 
  foo %>% 
  mutate(DVN = DV/DOSE) %>%
  group_by(IREP,TIME) %>% 
  summarise(med = median(DVN), lo = quantile(DVN,0.05), hi = quantile(DVN,0.95)) %>% 
  group_by(TIME) %>% summarise_at(vars(lo,med,hi),funs(median)) %>% ungroup

##' 
##' ## The mrgsolve data
##' 
sum_mrg <- 
  out %>%
  mutate(DVN = DV/DOSE) %>%
  group_by(IREP,TIME) %>%
  summarise(med = median(DVN), lo = quantile(DVN,0.05), hi = quantile(DVN,0.95)) %>%
  group_by(TIME) %>% summarise_at(vars(lo,med,hi),funs(median)) %>% ungroup


##' 
##' # Plot 1
##' 
##' Were are plotting the median 5th/50th/95th percentiles of the 
##' simulated data for both mrgsolve and nonmem
##' 
mrg <- sum_mrg %>% tidyr::gather(variable,value,lo:hi) %>% mutate(tool = "mrg")
non <- sum_nm %>%  tidyr::gather(variable,value,lo:hi) %>% mutate(tool = "nonmem")

#+ 
ggplot() + ggtitle("Lines: mrgsolve, Points: nonmem") + 
  geom_line(data = mrg, aes(TIME, value, col = variable, group = variable), lwd = 1) +
  geom_point(data = non, aes(TIME,value),col = "black", size = 2) + 
  scale_color_brewer(palette = "Set2", labels = c("95th", "5th", "50th")) 


##' # Plot 2
##' 
##' Look at 95% CI around the 5th/median/95th percentiles
##' 
##' 
#+
nm_q <- 
  foo %>% 
  mutate(DVN = DV/DOSE) %>%
  group_by(IREP,TIME) %>%  
  summarise(med = median(DVN), lo = quantile(DVN,0.05), hi = quantile(DVN,0.95)) %>%
  ungroup() %>%
  gather(metric,value,med:hi) %>% 
  group_by(metric,TIME) %>% 
  summarise(med = median(value), lo = quantile(value,0.025), hi = quantile(value,0.975)) %>%
  filter(TIME > 0)

#+
mrg_q <- 
  out %>%
  mutate(DVN = DV/DOSE) %>%
  group_by(IREP,TIME) %>%
  summarise(med = median(DVN), lo = quantile(DVN,0.05), hi = quantile(DVN,0.95)) %>%
  ungroup() %>%
  gather(metric,value,med:hi) %>% 
  group_by(metric,TIME) %>% 
  summarise(med = median(value), lo = quantile(value,0.025), hi = quantile(value,0.975)) %>%
  filter(TIME > 0)

mrg_q2 <- gather(mrg_q, q, value, med:hi)
nm_q2 <-  gather(nm_q,  q, value, med:hi)

##' Comparing 95% CI around 5th/50th/95th percentiles from
##' mrgsolve (ribbon and points) and nonmem (lines)
ggplot() + 
  geom_point(data = nm_q, aes(TIME, med, group = metric), size = 2) +
  geom_line(data  = mrg_q2, aes(TIME, value, group = paste(metric,q)), lwd = 0.8) +
  geom_ribbon(data=nm_q, aes(TIME, ymin = lo, ymax = hi, group = metric, fill = metric), alpha = 0.5) +
  scale_y_continuous(trans = "log10", breaks = 10^seq(-4,4)) + 
  ylab("DV") + ggtitle("Shaded area & points: 95% CI, nonmem; Solid lines: 95% CI, mrgsolve") + 
  scale_fill_brewer(palette = "Set2", labels = c("95th", "5th" ,"50th")) + 
  theme(legend.position = "top")



##' # Models
##' 
##' ## NONMEM  (104.ctl)
##' 
cat(readLines("model/104.ctl"), sep="\n")

##' ## mrgsolve
cat(mod@code, sep = "\n")

##' # Environment
devtools::session_info()
