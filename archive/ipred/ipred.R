
#' ---
#' title: IPRED
#' author: ""
#' date: ""
#' output: github_document
#' ---


.libPaths("/data/Rlibs")
library(mrgsolve)
library(tidyverse)


#' # Simulate the data with mrgsolve
#' 
#' There are no random effects in this model (RUV, IIV)
#' 

mod <- modlib("pk1") %>% param(CL = 1.1, V = 20, KA = 1.5) 
out <- mrgsim(mod, events = ev(amt = 100), carry_out = "amt,ii,addl,evid,cmt", delta = 0.5, recsort = 3)


plot(out, CP ~ time)

#' # Write out a data set for NONMEM

nm <- 
  out %>% 
  as_tibble() %>%
  mutate(C = '.') %>% 
  select(C,everything(), -EV, -CENT, -CP) %>% 
  mutate(DV  = '.', cmt = ifelse(evid==0,2,1))


head(nm)

write_csv(nm, path = "1001.csv")


#' # Simulate the data with NONMEM
#' 
#' Like the mrgsolve setup, no random effects (RUV or IIV)
#' 
metrumrg::NONR(1001, project = ".", 
               command = "/opt/NONMEM/nm74/nmqual/autolog.pl", 
               checkrunno=FALSE)

tab <- read_table("1001/TAB", skip = 1) %>% as.data.frame

#' ## Check
#' IPRED / PRED / A2/V are all the same
#' 
head(out)

head(tab)

summary(out$CP-tab[["IPRED"]])

#' # Simulate from NONMEM with IIV
#' 
#' Just adding IIV
#' 
metrumrg::NONR(1002, project = ".", 
               command = "/opt/NONMEM/nm74/nmqual/autolog.pl", 
               checkrunno=FALSE)

tab2 <- read_table("1002/TAB", skip = 1)

head(tab2)


#' # Compare IIV with no IIV
#' 
head(tab)

head(tab2)
