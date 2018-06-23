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

to_data_set <- function(x) {
  x <- as.data.frame(x)
  x <- mutate(x, C = '.', DV = '.', cmt = if_else(cmt==0, 2, cmt))
  dplyr::select(x, "C", everything())
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
sim <- function(x,e,...) {
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

##' # BOLUS
e <- ev(amt = 100)
out <- sim(mod, e)
data <- to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)


##' ## SUMMARY
summary(out$CP-outt$CP)


##' # INFUSION
e <- ev(amt = 100, rate = 10)
out <- sim(mod, e)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # INFUSION, LAG
e <- ev(amt = 100, rate = 10, cmt = 2, LAGT = 5)
out <- sim(mod, e)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, MULTIPLE
e <- ev(amt = 100, rate = 10, cmt = 2, ii = 8, addl = 2)
out <- sim(mod, e)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # INFUSION, MULTIPLE, LAG, II < DUR
e <- ev(amt = 100, rate = 100/10, cmt = 2, ii = 6, addl = 2, LAGT = 2)
out <- as.data.frame(sim(mod, e, end = 48))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)


##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, MULTIPLE, LAG, II > DUR
e <- ev(amt = 100, rate = 10, cmt = 2, ii = 12, addl = 2, LAGT = 2)
out <- as.data.frame(sim(mod, e, end = 48))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # INFUSION, MULTIPLE, II > DUR
e <- ev(amt = 100, rate = 100/10, cmt = 2, ii = 24, addl = 2)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, MULTIPLE, II < DUR
e <- ev(amt = 100, rate = 100/10, cmt = 2, ii = 6, addl = 4)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)



##' # INFUSION, BIOAV
e <- ev(amt = 100, rate = 100/5, cmt = 2, ii = 12, addl = 4, BIOAV = 0.7)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, BIOAV, DUR
e <- ev(amt = 100, rate = -2, cmt = 2, ii = 12, addl = 4, BIOAV = 0.7, DUR2 = 7, MODE = 2)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)

# start <- seq(0,4*12,12)
# end <- start + 7
# summary(out$CP - outt$CP)
# ggplot(out, aes(time,CP)) + geom_line() + 
#   geom_vline(xintercept = c(start,end), col="firebrick")


##' # INFUSION, BIOAV, RATE
e <- ev(amt = 100, rate = -1, cmt = 2, ii = 12, addl = 4, BIOAV = 0.7, R2 = 100/9, MODE = 1)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, BIOAV, RATE, LAGTIME
e <- ev(amt = 100, rate = -1, cmt = 2, ii = 12, addl = 4,  RAT2 = signif(100/9,5), MODE = 1, LAGT = 3)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)

##' ## SUMMARY
summary(out$CP - outt$CP)

#+
start <- 3+seq(0,4*12,12)
end <- (start + 9)
summary(out$CP - outt$CP)
ggplot(out, aes(time,CP)) + geom_line() + 
  geom_vline(xintercept = c(start,end), col="firebrick")


##' # INFUSION, SS
e <- ev(amt = 100, rate = -1, cmt = 2, ii = 12, addl = 4,  RAT2 = signif(100/9,5), MODE = 1, ss = 1)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
out$NM <- outt$CP
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, SS, BIOAV
e <- ev(amt = 100, rate = -1, cmt = 2, ii = 12, addl = 4,  RAT2 = 100/10, MODE = 1, ss = 1, BIOAV = 0.3)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
out$NM <- outt$CP
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, SS, BIOAV, DUR
e <- ev(amt = 100, rate = -2, cmt = 2, ii = 12, addl = 4,  DUR2 = 2, MODE = 2, ss = 1, BIOAV = 0.3)
out <- sim(mod, e, end = 72)
plot(out)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # INFUSION, SS, BIOAV, DUR, DUR == II
e <- ev(amt = 100, rate = -2, cmt = 2, ii = 6, addl = 4,  DUR2 = 6, MODE = 2, ss = 1, BIOAV = 0.4)
out <- sim(mod, e, end = 72)
plot(out)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)

##' # INFUSION, SS, BIOAV, DUR, II multiple of DUR
e <- ev(amt = 100, rate = -2, cmt = 2, ii = 6, addl = 4,  DUR2 = 12, MODE = 2, ss = 1, BIOAV = 0.4)
out <- sim(mod, e, end = 72)
plot(out)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # INFUSION, SS, BIOAV, DUR,  II < DUR
e <- ev(amt = 100, rate = -2, cmt = 2, ii = 6, addl = 4,  DUR2 = 10, MODE = 2, ss = 1, BIOAV = 0.65)
out <- sim(mod, e, end = 72)
plot(out)
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # BOLUS, SS
e <- ev(amt = 100, cmt = 2, ii = 12, addl = 4, ss = 1)
out <- as.data.frame(sim(mod, e, end = 72))
data = to_data_set(out)
sv(data, "data/101.csv")
outt <- run(101)
out$NM <- outt$CP
head(out)

##' ## SUMMARY
summary(out$CP - outt$CP)


##' # Control stream
#+ comment = "  "
writeLines(readLines("model/101.ctl"))
#+
devtools::session_info()

