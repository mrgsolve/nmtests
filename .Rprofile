
mpn_date <- "2022-03-21"

mpn_url <- paste0("https://mpn.metworx.com/snapshots/stable/", mpn_date, "/")

options(
  repos = c(
    MPN = mpn_url, 
    # a value must be set to CRAN or R will complain, so we'll point both to MPN
    CRAN = mpn_url
  )
)

r_version <- "4.1"

correct_r <- grepl(
  paste0("R version ", r_version),
  R.version[["version.string"]], 
  fixed = TRUE
)

# source after setting the repos to make sure renv will see those repo versions
if (correct_r) {
  source("renv/activate.R")
} else {
  stop(paste0("The project only works with R ", r_version))
}

rm(mpn_date, mpn_url, r_version, correct_r)
