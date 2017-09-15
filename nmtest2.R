##' ---
##' output: pdf_document
##' title: "INFUSION TESTS"
##' ---


Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio-server/bin/pandoc")


.libPaths("/data/Rlibs")
library(mrgsolve)
library(dplyr)
library(readr)
library(ggplot2)

carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")

to_data_set <- function(x,id = NULL) {
  x <- as.data.frame(x)
  x <- mutate(x, C = '.', DV = '.', cmt = if_else(cmt==0, 2, cmt))
  x <- dplyr::select(x, "C", everything())
  if(is.numeric(id)) x <- mutate(x,ID = id)
  x
}
sv <- function(x,file) {
  write.csv(file = file, row.names = FALSE, quote = FALSE, x)
}
nonr <- metrumrg::NONR
run <- function(number) {
  nonr(number, project = "model", command = "/opt/NONMEM/nm73/nmqual/autolog.pl", checkrunno=FALSE)
  return(tabread(number))
}
tabread <- function(number) {
  tab <- file.path("model", number, "TAB")
  if(file.exists(tab)) return(read_table(tab, skip=1))
  stop("the run failed")
}
sim <- function(x, e,...) {
  mrgsim(x, events = e, carry.out = carry, digits = 5, ...) 
}


code <- '
$SET req = ""
$PARAM CL = 1, V = 30, KA = 1.5
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
mod <- update(mod, end=72)


##' ### Bolus
ev <- ev(amt = 100, ii = 12, addl = 3) 
out1 <- sim(mod,ev)
data1 <- to_data_set(out1, 1)

##' ### Infusion
ev <- ev(amt = 100, ii = 12, addl =3, rate = 100/10, cmt = 2) 
out2 <- sim(mod,ev)
data2 <- to_data_set(out2, 2)

##' ### Infusion, lag
ev <- ev(amt = 100, ii = 12, addl=3, rate = 100/10, LAGT = 5, cmt = 2) 
out3 <- sim(mod,ev)
data3 <- to_data_set(out3, 3)

##' ### Infusion, lag, bioav
ev <- ev(amt = 100, ii = 12, addl=3, rate = 100/10, LAGT = 5, BIOAV = 0.412, cmt = 2) 
out4 <- sim(mod,ev)
data4 <- to_data_set(out4, 4)

##' ### Infusion,  bioav, ss
ev <- ev(amt = 100, ii = 12, addl=3, rate = 100/10, LAGT = 0, BIOAV = 0.412, ss = 1, cmt = 2) 
out5 <- sim(mod,ev)
data5 <- to_data_set(out5, 5)

##' ### Infusion, bioav, ss, II < DUR
ev <- ev(amt = 100, ii = 6, addl = 3, rate = 100/10, BIOAV = 0.812, ss = 1, cmt = 2) 
out6 <- sim(mod,ev)
data6 <- to_data_set(out6, 6)

##' ### Infusion, ss, II < DUR
ev <- ev(amt = 100, ii = 6, addl = 3, rate = 100/10, ss = 1, cmt = 2) 
out6.1 <- sim(mod,ev)
data6.1 <- to_data_set(out6.1, 6.1)

##' ### Infusion, ss, II multiple of DUR
ev <- ev(amt = 100, ii = 6, addl = 3, rate = signif(100/12,5), ss = 1, cmt = 2) 
out6.2 <- sim(mod,ev)
data6.2 <- to_data_set(out6.2, 6.2)

##' ### Infusion, bioav, ss, II == DUR
ev <- ev(amt = 100, ii = 10, addl=3, rate = 100/10, LAGT = 0, BIOAV = 0.412, ss = 1, cmt = 2) 
out7 <- sim(mod,ev)
data7 <- to_data_set(out7, 7)

##' ### Infusion,, ss, II == DUR
ev <- ev(amt = 100, ii = 10, addl=3, rate = 100/10, LAGT = 0, ss = 1, cmt = 2) 
out7.1 <- sim(mod,ev)
data7.1 <- to_data_set(out7.1, 7.1)

##' ### Bolus, bioav, ss
ev <- ev(amt = 100, ii = 12, addl=3,  LAGT = 0, BIOAV = 0.412, ss = 1) 
out8 <- sim(mod,ev)
data8 <- to_data_set(out8, 8)

##' ### Bolus, lag, bioav
ev <- ev(amt = 100, ii = 12, addl=3,  LAGT = 5, BIOAV = 0.412) 
out9 <- sim(mod,ev)
data9 <- to_data_set(out9, 9)

##' ### Infusion / bolus
ev <- ev(amt = 100, rate = 10) + ev(time = 12, amt = 50)
out10 <- sim(mod,ev)
data10 <- to_data_set(out10, 10)

##' ### Infusion (D_) lag, BIOAV
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 12, addl = 5, LAGT = 5, BIOAV = 0.61)
out11 <- sim(mod,ev)
data11 <- to_data_set(out11,11)

##' ### Infusion (D_) ss, BIOAV
ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 12, addl = 5, ss = 1, BIOAV = 0.61)
out12 <- sim(mod,ev)
data12 <- to_data_set(out12,12)



##' ## Simulate
sims <- list(out1,out2,out3,out4,out5,out6,out6.1,out6.2,out7,out7.1,
             out8,out9,out10,out11,out12)
sims <- lapply(sims, as.data.frame)
sims <- bind_rows(sims)
data <- bind_rows(data1,data2,data3,data4,data5,data6,data6.1,data6.2,data7,data7.1,
                  data8,data9,data10,data11,data12)
sv(data, "data/101.csv")
out <- run(101)


##' # Overall Summary
dim(out)
dim(sims)
summary(out$CP - sims$CP)


data$NM <- out$CP
data$MRGSIM <- sims$CP

##' # Summary by RUN
group_by(data,ID) %>% 
  mutate(diff = NM - MRGSIM) %>%
  summarise(mean = mean(diff), max = max(diff), min = min(diff))

##' # Plot
ggplot(data = data) + 
  geom_point(aes(time,NM),color = "firebrick") + 
  geom_line(aes(time,MRGSIM,group = ID)) +
  facet_wrap(~ID)


##' # Control stream
#+ comment = "  "
writeLines(readLines("model/101.ctl"))
#+
devtools::session_info()



