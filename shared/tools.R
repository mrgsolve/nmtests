#' # Functions
#' 
#' These functions assemble data sets, run simulations, and gather outputs. All 
#' scenarios are handled in exactly the same way.
#' 
#' ## Save `mrgsim` output as a `nonmem` input data set
to_data_set <- function(x, id = NULL) {
  x <- as.data.frame(x)
  x <- mutate(x, C = '.', DV = '.', cmt = if_else(cmt==0, 2, cmt))
  x <- select(x, "C", everything())
  if(is.numeric(id)) x <- mutate(x,ID = id)
  x
}

#' ## Save the `nonmem` input data set
fsave <- function(x,file) {
  fwrite(file = file, quote = FALSE, x = x, na = '.')
}

#' ## Run `nonmem`
psn_execute <- function(number, project = '.') {
  #execute -model_subdir -silent -directory=run1001 model/1001.ctl
  directory <- paste0("-directory=", number)
  target <- file.path(project, paste0(number, ".ctl"))
  args <- c("-silent", "-model_subdir", directory, target)
  system2("execute", args=args)
  return(tabread(number))
}

#' ## Read in `nonmem` simulation results
tabread <- function(number) {
  tab <- file.path(number, "TAB")
  if(file.exists(tab)) return(fread(tab, skip=1))
  stop("the run failed")
}

#' ## Simulate a scenario with `mrsim()`
sim <- function(x, e,...) {
  carry <- c("cmt", "amt","ii", "addl", "rate", "evid", "ss")
  mrgsim(x, events = e, carry_out = carry, digits = 5, recsort = 3, ...) 
}

push_back <- function(env, ev, descr) {
  n <- length(env$ev)+1
  m <- length(env$descr)+1
  env$ev[[n]] <- ev
  env$descr[[m]] <- descr
}

comp_plot <- function(comp) {
  id <- comp$ID[1]
  title <- paste0("ID: ", id, "; line: mrgsolve, point: NONMEM")
  ggplot(data = comp) + 
    geom_point(aes(time, NONMEM),color = "firebrick") + 
    geom_line(aes(time, MRGSIM, group = ID)) +
    ylab("Simulated value") + xlab("Time") + 
    ggtitle(label = NULL, subtitle = title) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(0, 130, 24))  
}
