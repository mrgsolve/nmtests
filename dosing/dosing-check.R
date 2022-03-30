mrgsolve.loc <- NULL
library(dplyr)
library(mrgsolve, lib.loc = mrgsolve.loc)
library(data.table)
library(testthat)
library(here)

setwd(here("dosing"))

key <- fread("results/1001-run-key.csv")
mod <- mread("1001.mod")
data <- fread("data/1001.csv", na.strings = '.')
tab <- fread("1001/TAB", skip = 1, na.strings='.')
out <- mrgsim_df(mod, data, recsort = 3, digits = 5)
out <- left_join(out, key, by = "ID")

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
