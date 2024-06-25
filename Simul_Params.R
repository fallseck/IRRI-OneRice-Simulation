# IRRI OneRice Breeding Strategies
# Assuming a complex Single trait represented by Grain Yield

# Simulation parameters -----
library(AlphaSimR)
scenarioName = c("Baseline","5-Years","3-Years","2-Years-WP","2-Years-BP") 
source("Suppl_Functions .R")

# Specie parameters ---------
nChr = 12
nQTL = 250
nSNP = 100
Ne = 60
genLen = 1.20
Bp = 3e+07 

# Burnin parameters ---------
CohortBrn = 6
nFounders = 80
nParentsBrn = 80
nCrossBrn = 100
nProgBrn = 120
nStg1 = 25
nStg2 = 600
nStg3 = 50
nBrn = 40#

# Breeding parameters -------
Cohort5Y = 5
Cohort3Y = 4
Cohort2Y = 3
iBreed = nBrn+1
nBreed = 30 
nYears = nBrn+nBreed
nParents = 40
nCross = 30
nProg = 240
nFamOYT = 40
nFamEST = 7
nAYT = 30 

# Variance components -------
genMean = 0 
genVar = 1
varGxY = c(0,4,8) #c(1,2,3)
plantH2 = 0.001
plotH2 = 0.10
Fixeff = 1L
