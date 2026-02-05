library(dplyr)
library(mrgsolve)
library(data.table)
library(here)


# C ID TIME EVID AMT CMT SS II ADDL RATE DV

data <- as_data_set(evd(amt = 100, ii = 24, addl = 2, rate = 0, ss =0))
data <- expand_observations(data, times = c(0.5,1,2,4,8,12,24))
data <- mutate(data, C = NA, DV = NA, ID = 1)
data <- select(data, C, ID, TIME, EVID, AMT, CMT, SS, II, ADDL, RATE, DV)

fwrite(data, file = here("dosing/data/8888.csv"), na = '.')


