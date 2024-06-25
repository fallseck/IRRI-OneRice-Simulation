m = Sys.getenv("SLURM_ARRAY_TASK_ID")
m = as.integer(m)

nReps=1

for (REP in 1:nReps) {
  cat("----------------------------------------------------\n")
  cat(paste0("ITERATION ", m, " [][][][][] \n"))
  
  # Create data.frame for results
  Mean = data.frame(Cycle = numeric(nBreed+1),
                    Rep = rep(m, c(nBreed+1)),
                    GxE = rep(g, c(nBreed+1)),
                    Scenario = numeric(nBreed+1),
                    meanG_Sg1 = numeric(nBreed+1),
                    varG_Sg1 = numeric(nBreed+1),
                    meanG_Sg2 = numeric(nBreed+1),
                    varG_Sg2 = numeric(nBreed+1),
                    Accu = numeric(nBreed+1),
                    Kind = numeric(nBreed+1),
                    stringsAsFactors=FALSE)
  
  Herr = data.frame(Year = 1:nYears,
                    REP = rep(m, nYears),
                    GxE = rep(g, nYears),
                    Scenario = numeric(nYears),
                    LST = numeric(nYears),
                    OYT = numeric(nYears),
                    EST = numeric(nYears),
                    AYT = numeric(nYears),
                    Kind = numeric(nYears),
                    stringsAsFactors=FALSE)
  
  tmp = data.frame(Year = 0:nYears,
                   REP = rep(m, nYears+1),
                   GxE = rep(g, nYears+1),
                   Scenario = numeric(nYears+1),
                   stringsAsFactors=FALSE)
  
  qtl = alphaP = matrix(NA, nrow = nYears+1, ncol = nQTL*nChr)
  colnames(qtl) = colnames(alphaP) = paste0('qtl_', 
                                                                 (1:(nQTL*nChr)))
  snp = matrix(NA, nrow = nYears+1, ncol = nSNP*nChr)
  colnames(snp) = paste0('snp_', (1:(nSNP*nChr)))
  param = tibble::tibble(Param = c("QTL", "SNP" ,"Ne","MeanG", "VarG","GxE",
                                   "Plant_H2","Plot_H2","nYears", "nBurnin"),
                         Value = c(nQTL, nSNP, Ne, genMean, genVar,g, 
                                   plantH2, plotH2, nYears, nBrn))
  
  cat("----------------------------------------------------\n")
  cat(paste0("SIMULATION OF FOUNDER HAPLOTYPES ... \n"))
  cat("----------------------------------------------------\n")
  
  # Simulate founder haplotypes
  founderPop = runMacs2(nInd = nFounders,
                        nChr = nChr,
                        segSites = nQTL+nSNP,
                        genLen = genLen,
                        Ne = Ne,
                        bp = Bp)
  # Set simulation parameters
  SP = SimParam$new(founderPop)
  
  # Add trait (assuming single trait as grain yield)
  SP$addTraitAG(nQtlPerChr = nQTL,
                mean = genMean,
                var = genVar,
                varGxE = g)
  
  # Add SNP Chip
  SP$addSnpChip(nSnpPerChr = nSNP)
  
  # Initiate the founders population
  FounderPop = newPop(founderPop)
  
  #Frequencies
  qtl[1,] = allFreq(pullQtlGeno(founderPop))$maf;alphaP[1,] = alphaFreq(pullQtlGeno(founderPop))$p
  snp[1,] = allFreq(pullSnpGeno(founderPop))$maf
  tmp$Scenario[1] = "BURNIN"; tmp$Kind[1] = "PS"
  
  # FILL THE BREEDING PROGRAM PIPELINE 
  
  cat(paste0("FILLING BREEDING PIPELINE ... \n"))
  cat("----------------------------------------------------\n")
  
  pVal = runif(CohortBrn)
  
  for (i in 1:CohortBrn) {
    
    # Crossing block and F1s Validation: 2 Seasons (Year 1)
    F1 = randCross(pop = FounderPop,
                   nCrosses = nCrossBrn,
                   nProgeny = 1)
    
    # RGA system: 2 Seasons (Year 2)
    if(i < 6){
      F2 = self(pop = F1,
                nProgeny = nProgBrn)
      F3 = self(F2)
      F4 = self(F3)
      F5 = self(F4)
      F6 = self(F5)
    }
    
    # Line testing and MAS: 1 Season (Year 3)
    if (i < 5){
      F6_LST = setPheno(pop = F6,
                        H2 = plantH2,
                        reps = 1, #1 Location
                        fixEff = Fixeff,
                        p = pVal[i+2])
    }
    
    # 1st Stage of Yield Evaluation: 1 Season (Year 4)
    if (i < 4){
      F7_Stg1 = selectWithinFam(pop = F6_LST,
                                nInd = nStg1,
                                use = 'pheno')
      F7_Stg1 = self(F7_Stg1)
      F7_Stg1 = setPheno(pop = F7_Stg1,
                         H2 = plantH2,
                         reps = 1, #1 Location
                         fixEff = Fixeff,
                         p = pVal[i+3])
    }
    
    # 2nd Stage of Yield Evaluation: 2 Seasons (Year 5)
    if (i < 3){
      F8_Stg2 = selectInd(pop = F7_Stg1,
                          nInd = nStg2,
                          use = 'pheno')
      F8_Stg2 = self(F8_Stg2)
      F8_Stg2 = setPheno(pop = F8_Stg2,
                         H2 = plotH2,
                         reps = 4, #4 Locations
                         fixEff = Fixeff,
                         p = pVal[i+4])
    }
    
    # 3rd Stage of Yield Evaluation: 2 Seasons (Year 6)
    if (i < 2){
      F9_Stg3 = selectInd(pop = F8_Stg2,
                          nInd = nStg3,
                          use = 'pheno')
      F9_Stg3 = self(F9_Stg3)
      F9_Stg3 = setPheno(pop = F9_Stg3,
                         H2 = plotH2,
                         reps = 4, #4 Locations
                         fixEff = Fixeff,
                         p = pVal[i+5])
    }
    
  }#END Pipeline
  
  TrainPop = newEmptyPop()
  P = runif(nYears)
  
  for (b in 1:nBrn) {
    ib = b+1
    # Select parents from previous evaluated populations
    GermplasmPool = c(F8_Stg2,F9_Stg3)
    
    ParentsBrn = selectInd(pop = GermplasmPool,
                           nInd = nParentsBrn,
                           use = 'pheno')
    
    # 3rd Stage of Yield Evaluation: 2 Seasons (Year 6)
    F9_Stg3 = selectInd(pop = F8_Stg2,
                        nInd = nStg3,
                        use = 'pheno')
    F9_Stg3 = self(F9_Stg3)
    F9_Stg3 = setPheno(pop = F9_Stg3,
                       H2 = plotH2,
                       #h2 = ploth2,
                       reps = 4, #4 Locations
                       fixEff = Fixeff,
                       p = P[b])
    
    # 2nd Stage of Yield Evaluation: 2 Seasons (Year 5)
    F8_Stg2 = selectInd(pop = F7_Stg1,
                        nInd = nStg2,
                        use = 'pheno')
    F8_Stg2 = self(F8_Stg2)
    F8_Stg2 = setPheno(pop = F8_Stg2,
                       H2 = plotH2,
                       reps = 4, #4 Locations
                       fixEff = Fixeff,
                       p = P[b])
    
    # 1st Stage of Yield Evaluation: 1 Season (Year 4)
    F7_Stg1 = selectWithinFam(pop = F6_LST,
                              nInd = nStg1,
                              use = 'pheno')
    F7_Stg1 = self(F7_Stg1)
    F7_Stg1 = setPheno(pop = F7_Stg1,
                       H2 = plantH2,
                       reps = 1, #1 Location
                       fixEff = Fixeff,
                       p = P[b])
    
    # Line testing and MAS: 1 Season (Year 3)
    F6_LST = setPheno(pop = F6,
                      H2 = plantH2,
                      reps = 1, #1 Location
                      fixEff = Fixeff,
                      p = P[b])
    
    # RGA system: 2 Seasons (Year 2)
    F2 = self(pop = F1,
              nProgeny = nProgBrn)
    F3 = self(F2)
    F4 = self(F3)
    F5 = self(F4)
    F6 = self(F5)
    
    # Crossing block and Hybridity validation: 2 Seasons (Year 1)
    F1 = randCross(pop = ParentsBrn,
                   nCrosses = nCrossBrn,
                   nProgeny = 1)
    
    # Build Training Pop for prediction from Yield trial over last 3 years
    if(b > (nBrn-3)){
      TrainPop = c(TrainPop,
                   selectInd(GermplasmPool,200)) 
    }
    
    cat(paste0("Saving Data From Crossing Year ", b,  " ....\n"))
    
    # Save results
   
    #LSTs
    Herr$LST[b] = varG(F6_LST)[1]/varP(F6_LST)[1]
    # Stage 1 (F7s)
    Herr$OYT[b] = varG(F7_Stg1)[1]/varP(F7_Stg1)[1]
    # Stage 2 (F8s)
    Herr$EST[b] = varG(F8_Stg2)[1]/varP(F8_Stg2)[1]
    # Stage 3 (F9s)
    Herr$AYT[b] = varG(F9_Stg3)[1]/varP(F9_Stg3)[1]
    
    #AllFreq
    qtl[ib,] = allFreq(pullQtlGeno(ParentsBrn))$maf
    alphaP[ib,] = alphaFreq(pullQtlGeno(ParentsBrn))$p
    snp[ib,] = allFreq(pullSnpGeno(ParentsBrn))$maf
    tmp$Scenario[ib] = "BURNIN"; tmp$Kind[ib] = "PS"
    
    if(b==nBrn){
      cat("----------------------------------------------------\n")
      cat(paste0("END OF THE BURN-IN PHASE \n"))
    }
  }#END of the BurnIn
  
  ECP = selectInd(GermplasmPool,80)
  #AllFreq
  qtl[ib,] = allFreq(pullQtlGeno(ECP))$maf;alphaP[ib,] = alphaFreq(pullQtlGeno(ECP))$p
  snp[ib,] = allFreq(pullSnpGeno(ECP))$maf
  Mean$Cycle[1] = 0; Mean$meanG_Sg1[1] = Mean$meanG_Sg2[1] = meanG(ECP)
  Mean$varG_Sg1[1] = Mean$varG_Sg2[1] = varG(ECP)[1]; Mean$Accu[1] = cor(gv(ECP),EBV(ECP))
  pVal = runif(Cohort5Y)
  save.image(paste0("BURNIN_OneRice_",g,"_",m,".RData"))
  
}#END of ITERATIONS
