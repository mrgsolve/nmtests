
-   [1 Input data](#1-input-data)
    -   [1.1 Constant dose](#11-constant-dose)
    -   [1.2 Time-varying dose](#12-time-varying-dose)
-   [2 mrgsolve - ode](#2-mrgsolve---ode)
-   [3 Constant dose](#3-constant-dose)
    -   [3.1 ODE](#31-ode)
    -   [3.2 PRED](#32-pred)
-   [4 Time-varying dose](#4-time-varying-dose)
    -   [4.1 ODE](#41-ode)
    -   [4.2 PRED](#42-pred)
-   [5 Plots](#5-plots)
    -   [5.1 Constant dose](#51-constant-dose)
    -   [5.2 Time-varying](#52-time-varying)
    -   [5.3 Time-varying, ODE,
        NM/mrgsolve](#53-time-varying-ode-nmmrgsolve)
-   [6 Spot check analytical
    solution](#6-spot-check-analytical-solution)

``` r
library(dplyr)
library(data.table)
library(here)
library(ggplot2)
theme_set(theme_bw())
library(mrgsolve)

source(here("shared/tools.R"))
setwd(here("claret"))
```

# 1 Input data

## 1.1 Constant dose

``` r
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
```

## 1.2 Time-varying dose

``` r
data2 <- mutate(
  data, 
  DOSE = ifelse(TIME == 0.5, 5, DOSE)
)

fwrite(data, file = here("claret/data/claret001.csv"), na = '.', quote = FALSE)
fwrite(data2, file = here("claret/data/claret002.csv"), na = '.', quote = FALSE)
```

# 2 mrgsolve - ode

``` r
mod <- mread(here("claret/ode.txt")) %>% zero_re()
```

    ## Building ode_txt ... done.

``` r
z <- mrgsim(mod, data2, output = "df")
```

# 3 Constant dose

## 3.1 ODE

``` r
a <- psn_execute("c001")
```

## 3.2 PRED

``` r
b <- psn_execute("c002")
```

# 4 Time-varying dose

## 4.1 ODE

``` r
c <- psn_execute("c003")
```

## 4.2 PRED

``` r
d <- psn_execute("c004")
```

# 5 Plots

## 5.1 Constant dose

``` r
ggplot() + 
  geom_point(data = a, aes(TIME, Y)) + 
  geom_line(data = b,  aes(TIME, Y), color = "firebrick") + 
  ggtitle("Point: ODE, Line: PRED")
```

![](claret_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## 5.2 Time-varying

``` r
ggplot() + 
  geom_point(data = c, aes(TIME, Y)) + 
  geom_line(data = d,  aes(TIME, Y), color = "firebrick") + 
  ggtitle("Point ODE; Line: PRED")
```

![](claret_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## 5.3 Time-varying, ODE, NM/mrgsolve

``` r
ggplot() + 
  geom_point(data = c, aes(TIME, Y)) + 
  geom_line(data = z,  aes(TIME, RESP), color = "firebrick") + 
  ggtitle("Point nonmem; Line: mrgsolve")
```

![](claret_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# 6 Spot check analytical solution

``` r
KG = 0.6
KS0 = 0.4
GAMMA = 0.8
BASE = 70
RESP_0 = BASE
DOSE = 5

TIME = 0.5
KS = KS0 *exp( -GAMMA * TIME);
RESP_0 * exp((KG*TIME) - (KS0 * log(DOSE)/GAMMA) * (1-exp(-GAMMA*TIME)))
```

    ## [1] 72.47166

``` r
filter(d, TIME==0.5)
```

    ##    ID TIME      Y
    ## 1:  1  0.5 72.472

``` r
TIME = 0.4
DOSE = 10
KS = KS0 *exp( -GAMMA * TIME);
RESP_0 * exp((KG*TIME) - (KS0 * log(DOSE)/GAMMA) * (1-exp(-GAMMA*TIME)))
```

    ## [1] 64.9237

``` r
filter(d, TIME==0.4)
```

    ##    ID TIME      Y
    ## 1:  1  0.4 64.924

``` r
TIME = 1.0
DOSE = 10
KS = KS0 *exp( -GAMMA * TIME);
RESP_0 * exp((KG*TIME) - (KS0 * log(DOSE)/GAMMA) * (1-exp(-GAMMA*TIME)))
```

    ## [1] 67.66112

``` r
filter(d, TIME==1.0)
```

    ##    ID TIME      Y
    ## 1:  1    1 67.661
