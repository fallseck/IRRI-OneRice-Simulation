rm(list = ls())
time1 = Sys.time()
source("Simul_Params.R")

for (g in varGxY) {
  cat("----------------------------------------------------\n")
  cat(paste0("[][][][][] GxE LEVEL ", g, " [][][][][] \n"))
  
  source("Burnin_IRRI_OneRice.R")
  source("Baseline_OneRice_IRRI.R")
  source("5Y_OneRice_IRRI.R")
  source("3Y_OneRice_IRRI.R")
  source("2Y_OneRice_IRRI.R")
  source("2YM_OneRice_IRRI.R")
}

time2 = Sys.time()

cat("----------------------------------------------------\n")
cat(paste0("END OF THE BREEDING PROGRAM \n"))
cat(paste0("SIMULATION RUNTIME: ", runTime(time1,time2), "\n"))
cat("----------------------------------------------------\n")