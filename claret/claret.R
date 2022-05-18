#' ---
#' output: 
#'   github_document: 
#'     toc: true
#'     number_sections: true
#' ---

#+ message = FALSE, warning = FALSE
library(dplyr)
library(data.table)
library(here)
library(ggplot2)

source(here("shared/tools.R"))
setwd(here("claret"))

#' 
#' # Input data
#' 
#' ## Constant dose
#' 
data <- expand.grid(
  C = NA_character_,
  ID = 1, 
  TIME = seq(0,1, 0.1), 
  DOSE = 10, 
  CMT = 1, 
  EVID=0, 
  MDV = 1, 
  DV = NA_real_
)

#' 
#' ## Time-varying dose
#' 
data2 <- mutate(
  data, 
  DOSE = ifelse(TIME == 0.5, 5, DOSE)
)

fwrite(data, file = here("claret/data/claret001.csv"), na = '.', quote = FALSE)
fwrite(data2, file = here("claret/data/claret002.csv"), na = '.', quote = FALSE)

#' # Constant dose
#' 
#' ## ODE
a <- psn_execute("c001")
#' ## PRED
b <- psn_execute("c002")

#' # Time-varying dose
#' 
#' ## ODE
c <- psn_execute("c003")
#' ## PRED
d <- psn_execute("c004")



#' # Plots
#' 
#' ## Constant dose
ggplot() + 
  geom_point(data = a, aes(TIME, Y)) + 
  geom_line(data = b,  aes(TIME, Y), color = "firebrick") + 
  ggtitle("Point: ODE, Line: PRED")

#' ## Time-varying
ggplot() + 
  geom_point(data = c, aes(TIME, Y)) + 
  geom_line(data = d,  aes(TIME, Y), color = "firebrick") + 
  ggtitle("Point ODE; Line: PRED")


#' 
#' # Spot check analytical solution
#' 
KG = 0.6
KS0 = 0.4
GAMMA = 0.8
BASE = 70
RESP_0 = BASE
DOSE = 5

TIME = 0.5
KS = KS0 *exp( -GAMMA * TIME);
RESP_0 * exp((KG*TIME) - (KS0 * log(DOSE)/GAMMA) * (1-exp(-GAMMA*TIME)))
filter(d, TIME==0.5)

TIME = 0.4
DOSE = 10
KS = KS0 *exp( -GAMMA * TIME);
RESP_0 * exp((KG*TIME) - (KS0 * log(DOSE)/GAMMA) * (1-exp(-GAMMA*TIME)))
filter(d, TIME==0.4)

TIME = 1.0
DOSE = 10
KS = KS0 *exp( -GAMMA * TIME);
RESP_0 * exp((KG*TIME) - (KS0 * log(DOSE)/GAMMA) * (1-exp(-GAMMA*TIME)))
filter(d, TIME==1.0)
