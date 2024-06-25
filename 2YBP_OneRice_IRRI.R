m = Sys.getenv("SLURM_ARRAY_TASK_ID")
m = as.integer(m)

nReps=1

for (REP in 1:nReps) {
  
  load(paste0("BURNIN_OneRice_",g,"_",m,".RData"))

  cat("----------------------------------------------------\n")
  cat(paste0("RUNING OF THE BREEDING SCHEME ", scenarioName[5]," ... \n"))
  cat("----------------------------------------------------\n")
  
  for (j in 1:Cohort2Y) {
    
    Fixeff = j
    
    #Select parents from Base population
    Parents = selectInd(pop = ECP,
                        nInd = nParents,
                        use = "rand")
    
    # Crossing Block & F1s Validation: Season 1 & 2 (Year 1)
    Cross = randCross(pop = Parents,
                      nCrosses = nCross,
                      nProgeny = 1)
    
    if (j < 3){ 
      # Line Fixation in RGA and MAS: Season 3 (Year 2)
      F2 = self(pop = Cross,
                nProgeny = nProg)
      F3 = self(F2)
      F4 = self(F3)
      F4 = setPheno(pop = F4,
                    H2 = plantH2,
                    reps = 1, #1 Location
                    fixEff = Fixeff, 
                    p = pVal[j+1])
    }
    
    if (j < 2){
      # Genomic Prediction and Seed Amplification: Season 4 (Year 2)
      F5 = selectWithinFam(pop = F4,
                           nInd = 40,
                           use = 'pheno') 
      F5 = self(F5)
      # Genomic Prediction
      GS = RRBLUP(pop = TrainPop)
      
      F5_EBV = setEBV(pop = F5,
                      solution = GS,
                      value = 'bv')
      F5 = setPheno(pop = F5,
                    H2 = plantH2,
                    reps = 1, #1 Location
                    fixEff = Fixeff, 
                    p = pVal[j+1])
    }
  }#END Pipeline
  
  Mean$Scenario = scenarioName[5]
  Mean$Kind = "GS"
  
  # RUN 30 Years of Target Breeding Scheme
  
  for (b in iBreed:nYears) {
    
    Fixeff = b
    
    iCycle = b-nBrn+1
    
    # Select new parents from last estimated breeding values
    newParents = selectInd(pop = F5_EBV,
                           nInd = nParents,
                           use = 'ebv')
    
    # 1st Stage of Yield Evaluation on EST: Season 5 & 6 (Year 3)
    F6_EST = selectWithinFam(pop = F5,
                             nInd = nFamEST,
                             use = 'rand')
    F6_EST = self(F6_EST)
    F6_EST = setPheno(pop = F6_EST,
                      H2 = plotH2,
                      reps = 4, #4 Locations
                      fixEff = Fixeff, 
                      p = P[b])
    
    # Update the training population by including the new Estimation set
    TrainPop = c(TrainPop, F6_EST)
    
    # 2nd Stage of Yield Evaluation on AYT: Season 7 & 8 (Year 4)
    F7_AYT = selectInd(pop = F5_EBV,
                       nInd = nAYT,
                       use = 'ebv')
    # Advancement
    F7_AYT = self(F7_AYT)
    F7_AYT = setPheno(pop = F7_AYT,
                      H2 = plotH2,
                      reps = 4, #4 Locations
                      fixEff = Fixeff, 
                      p = P[b])
    
    # Genomic Prediction and Seed Amplification: Season 4 (Year 2)
    F5 = selectWithinFam(pop = F4,
                         nInd = 40,
                         use = 'pheno') 
    F5 = self(F5)
    
    GS = RRBLUP(pop = TrainPop)
    
    F5_EBV = setEBV(pop = F5,
                    solution = GS,
                    value = 'bv')
    
    F5 = setPheno(pop = F5,
                  H2 = plantH2,
                  reps = 1, #1 Location
                  fixEff = Fixeff, 
                  p = P[b])
    
    # Line Fixation in RGA: Season 3 (Year 2)
    F2 = self(pop = Cross,
              nProgeny = nProg)
    F3 = self(F2)
    F4 = self(F3)
    F4 = setPheno(pop = F4,
                  H2 = plantH2,
                  reps = 1, #1 Location
                  fixEff = Fixeff, 
                  p = P[b])
    
    # Crossing Block & F1s Validation: Season 1 & 2 (Year 1)
    Cross = randCross(pop = newParents,
                      nCrosses = nCross,
                      nProgeny = 1)
    
    # Save results
    cat(paste0("Saving Data From Crossing Year ", b,  " ....\n"))
    
    # Stage 1 
    Mean$meanG_Sg1[iCycle] =  meanG(F5_EBV)
    Mean$varG_Sg1[iCycle] = varG(F5_EBV)[1]
    Mean$Accu[iCycle] = cor(gv(F5_EBV),ebv(F5_EBV))
    Mean$Cycle[iCycle] = iCycle-1
    # Stage 1 (ESTs)
    Herr$EST[b] = varG(F6_EST)[1]/varP(F6_EST)[1]
    
    #Stage 2 (AYTs)
    Mean$meanG_Sg2[iCycle] =  meanG(F7_AYT)
    Mean$varG_Sg2[iCycle] = varG(F7_AYT)[1]
    Herr$AYT[b] = varG(F7_AYT)[1]/varP(F7_AYT)[1]
    
    #AllFreq
    qtl[b+1,] = allFreq(pullQtlGeno(newParents))$maf
    alphaP[b+1,] = alphaFreq(pullQtlGeno(newParents))$p
    snp[b+1,] = allFreq(pullSnpGeno(newParents))$maf
  
    tmp$Scenario[b+1] = scenarioName[5]
    tmp$Kind[b+1] = "GS"
    
  }#END of Breeding Program
  
  freqQtltmp = cbind(tmp,qtl) 
  freqAlphaP = cbind(tmp,alphaP)
  freqSnptmp = cbind(tmp,snp)  
  
  saveRDS(list(PARAM=param, 
               MEAN=Mean, 
               H2=Herr,
               QTL=freqQtltmp,
               alphaP=freqAlphaP,
               SNP=freqSnptmp),
          paste0("RES_",scenarioName[5], "_IRRI-OneRice_",g,"_",m,".rds"))
  
}#END OF ITERATIONS
 