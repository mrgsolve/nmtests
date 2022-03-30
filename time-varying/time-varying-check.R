mrgsolve.loc <- NULL
library(dplyr)
library(mrgsolve, lib.loc = mrgsolve.loc)
library(data.table)
library(testthat)
library(here)

setwd(here("time-varying"))

mod <- mread("2001.mod")
data <- fread("data/2001.csv", na.strings = '.')
tab <- fread("2001/TAB", skip = 1, na.strings='.')
out <- mrgsim_df(mod, data, recsort = 3, digits = 5)
out$descr <- "time-varying weight"

summary(out$CP - tab$CP)

out <- mutate(out, NM = tab$CP)

sp <- split(out, out$ID)

out <- lapply(sp,  function(x) {
  title <- x$descr[[1]]
  test <- x$CP
  ref <- x$NM
  test_that(title, {
    expect_equal(test, ref, tolerance=1e-5)  
  })
})
