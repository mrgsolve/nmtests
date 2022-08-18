#' # Assemble the scenarios
#' 
#' There is a lot of code here.  See the [results](#results) section to see 
#' input data objects next to simulated data output from mrgsolve and NONMEM.
#' 
#' - Doses into `cmt` 2 are intravascular and doses into `cmt` 1 are 
#'   extravascular
#' - `LAGT` sets the dosing lag time
#' - `BIOAV` sets the bioavailability fraction
#' 
data_1001 <- function() {
  env <- new.env()
  env$ev <- list()
  env$descr <- list()
  
  #+ 
  ev <- ev(amt = 100, ii = 24, addl = 3) 
  push_back(env,ev, "Bolus with additional")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3, LAGT = 12.13, BIOAV = 2.23, cmt = 2) 
  push_back(env, ev,"Bolus with lag time and bioav")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, cmt = 2) 
  push_back(env,ev,"Infusion with additional")
  #+ 
  ev <- ev(amt = 480, ii = 0, addl = 0, rate = 10, cmt = 2, BIOAV = 0.5) 
  push_back(env,ev,"Infusion with bioav factor")
  #+ 
  ev <- ev(amt = 480, ii = 0, addl = 0, rate = -2, cmt = 2, BIOAV = 0.5, MODE = 2, DUR2 = 48) 
  push_back(env,ev,"Infusion with bioav factor and dur")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/12, cmt = 1) 
  push_back(env,ev,"Infusion doses to depot, with additional")
  #+
  ev <- ev(amt = 100, ii = 24, addl=3, rate = 100/10, LAGT = 4.15, cmt = 2) 
  push_back(env,ev,"Infusion doses, with additional and lag time")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.25, BIOAV = 0.412, cmt = 2) 
  push_back(env,ev,"Infusion doses, with lag time and bioav factor")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3, rate = 100/10, LAGT = 3.16, BIOAV = 0.412, ss = 1, cmt = 2) 
  push_back(env,ev,"Infusion doses at steady-state, with lag time and bioav factor")
  #+
  ev <- ev(amt = 100, ii = 12, addl = 4, rate = 100/50, BIOAV = 0.812, ss = 1, cmt = 2) 
  push_back(env,ev,"Infusion doses at steady-state, with bioav factor")
  #+
  ev <- ev(amt = 100, ii = 12, addl = 3, rate = 100/50, ss = 1, cmt = 2) 
  push_back(env,ev,"Infusion doses at steady state, with lag time")
  #+
  ev <- ev(amt = 100, ii = 6, addl = 12, rate = signif(100/12,5), ss = 1, cmt = 2) 
  push_back(env,ev,"Infusion doses at steady state, II < DUR, no bioav factor")
  #+
  ev <- ev(amt = 100, ii = 10, addl = 8, rate = 0.412*100/10,  BIOAV = 0.412, ss = 1, cmt = 2) 
  push_back(env,ev,"Infusion doses at steady state where II == DUR, with bioav factor")
  #+
  ev <- ev(amt = 100, ii = 10, addl = 8, rate = 100/10, ss = 1, cmt = 2) 
  push_back(env,ev,"Infusion doses at steady state, where II == DUR")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3,  LAGT = 4, BIOAV = 0.412, ss = 1, cmt = 2) 
  push_back(env,ev,"Bolus doses at steady state, with bioav factor and lag time")
  #+
  ev <- ev(amt = 100, ii = 24, addl = 3,  LAGT = 5, BIOAV = 0.412, cmt = 2) 
  push_back(env,ev,"Bolus doses with lag time and bioav factor")
  #+
  ev <- 
    ev(amt = 100, cmt = 2, LAGT = 1) + 
    ev(time = 13, amt = 50, ii = 24, addl = 2, rate = 24)
  push_back(env,ev,"Bolus then infusion")
  #+
  ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, LAGT = 5, BIOAV = 0.61)
  push_back(env,ev,"Infusion with modeled duration, lag time, and bioav factor")
  #+
  ev <- ev(amt = 100, rate = -2, DUR2 = 9, MODE = 2, cmt = 2, ii = 24, addl = 3, ss = 1, BIOAV = 0.61)
  push_back(env,ev,"Infusion with modeled duration, at steady state with bioav factor")
  #+
  ev <- 
    ev(amt = 100, ii = 12, addl = 2, rate = 50, BIOAV = 0.61) + 
    ev(amt = 120, evid = 4, time = 50, BIOAV = 0.5, ii = 12, addl = 3)
  push_back(env,ev,"Reset and dose (EVID 4) with additional")
  #+
  ev <- ev(amt = 100, ii = 12, evid = 4, ss = 1) 
  push_back(env,ev,"Reset and dose (EVID 4) with SS=1")
  #+
  ev <- 
    ev(amt = 100, ii = 12, addl = 3, rate = 50, BIOAV = 0.61) + 
    ev(amt = 0, evid = 3, time = 50, cmt = 2, BIOAV = 1) + 
    ev(amt = 120, ii = 16, addl = 2, time = 54, BIOAV = 1)
  push_back(env,ev,"Reset (EVID 3) with additional")
  
  #+
  ev <- 
    ev(amt = 100, ii = 12, addl = 3, rate = 50, BIOAV = 0.61, LAGT = 10) + 
    ev(amt = 0, evid = 3, time = 50, cmt = 2, BIOAV = 1, LAGT = 10) + 
    ev(amt = 120, ii = 16, addl = 2, time = 54, BIOAV = 1, LAGT = 10)
  push_back(env,ev,"Reset (EVID 3) with additional and lag")
  
  #+
  ev <- 
    ev(amt = 100, ii = 24, addl = 3, ss = 1)  + 
    ev(amt = 50,  ii = 24, addl = 3, ss = 2, time = 12)
  push_back(env,ev,"Steady state 1 and 2")
  #+ 
  ev <- ev(amt = 0, rate = 100, ss = 1)
  push_back(env,ev,"Steady state infusion")
  
  #+
  update_id <- function(ev,id) mutate(ev, ID = id)
  
  #+
  runs <- tibble(ev = env$ev, descr = env$descr)
  runs <- mutate(runs, ID = seq(n()))
  runs <- mutate(runs, ev = map2(ev, ID, update_id))
  
  select(runs, ID, descr, ev, everything())
}
