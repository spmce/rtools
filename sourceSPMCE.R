#  INSTALL DEPENDENCIES (IF NECCESSARY)

mp <- c("ggplot2", "circumplex", "MplusAutomation")
np <- mp[!(mp %in% installed.packages()[,"Package"])]
if(length(np)) install.packages(np)

rm(mp, np)

# SOURCE SPCME FUNCTIONS

source("https://raw.githubusercontent.com/spmce/rtools/master/runSPMCE.R")
source("https://raw.githubusercontent.com/spmce/rtools/master/readSPMCE.R")
source("https://raw.githubusercontent.com/spmce/rtools/master/plotSPMCE.R")
source("https://raw.githubusercontent.com/spmce/rtools/master/svalSPMCE.R")

message("\nLoading depencies:")

suppressMessages(library(MplusAutomation))
suppressMessages(library(ggplot2))
suppressMessages(library(circumplex))

