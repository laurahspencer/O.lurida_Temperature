list.of.packages <- c("ggplot2", "car", "dplyr", "reshape2", "plotly", "tidyr", "ggthemes", "lme4", "multcomp", "rcompanion", "gridExtra", "cowplot", "survival") #add new libraries here 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load all libraries 
lapply(list.of.packages, FUN = function(X) {
  do.call("require", list(X)) 
})

sessionInfo()

