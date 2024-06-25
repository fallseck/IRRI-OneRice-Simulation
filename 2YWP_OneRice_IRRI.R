m = Sys.getenv("SLURM_ARRAY_TASK_ID")
m = as.integer(m)

nReps=1

for (REP in 1:nReps) {
  
  load(paste0("BURNIN_OneRice_",g,"_",m,".RData"))
  
  cat("----------------------------------------------------\n")
  cat(paste0("RUNING OF THE BREEDING SCHEME ", scenarioName[4]," ... \n"))
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
      # Line Fixation in RGA: Season 3 (Year 2)
      F2 = self(pop = Cross,
                nProgeny = nProg)
      F3 = self(F2)
      F3 = setPheno(pop = F3,
                    H2 = plantH2,
                    reps = 1, #1 Location
                    fixEff = Fixeff, 
                    p = pVal[j+1])
    }
    
    if (j < 2){
      # 1st Stage of Yield Evaluation (MAS): Season 3 (Year 2)
      F4 = selectWithinFam(pop = F3,
                           nInd = 40,
                           use = 'pheno') 
      F4 = self(F4)
      
      # 1st Stage of Yield Evaluation (Estimation Set): Season 4 (Year 2)
      F5_EST = selectWithinFam(pop = F4,
                               nInd = nFamEST,
                               use = 'rand')
      F5_EST = self(F5_EST)
      F5_EST = setPheno(pop = F5_EST,
                        H2 = plotH2,
                        reps = 4, #4 Locations
                        fixEff = Fixeff,
                        p = pVal[j+1])
      # Genomic Prediction
      GS = RRBLUP(pop = F5_EST)
      
      F4_EBV = setEBV(pop = F4,
                      solution = GS,
                      value = 'bv')
      
    }
  }#END Pipeline
  
  Mean$Scenario = scenarioName[4]
  Mean$Kind = "GS"
  
  # RUN 30 Years of Target Breeding Scheme
  
  for (b in iBreed:nYears) {
    
    Fixeff = b
    
    iCycle = b-nBrn+1
    
    # Select new parents from last estimated breeding values
    newParents = selectInd(pop = F4_EBV,
                           nInd = nParents,
                           use = 'ebv')
    
    # 2nd Stage of Yield Evaluation on AYT: Season 5 & 6 (Year 3)
    F6_AYT = selectInd(pop = F4_EBV,
                       nInd = nAYT,
                       use = 'ebv')
    # Advancement
    F6_AYT = self(F6_AYT)
    F6_AYT = setPheno(pop = F6_AYT,
                      H2 = plotH2,
                      reps = 4, #4 Locations
                      fixEff = Fixeff, 
                      p = P[b])
    
    # 1st Stage of Yield Evaluation (MAS): Season 3 (Year 2)
    F4 = selectWithinFam(pop = F3,
                         nInd = 40,
                         use = 'pheno') 
    F4 = self(F4)
    
    # 1st Stage of Yield Evaluation (Estimation Set): Season 4 (Year 2)
    F5_EST = selectWithinFam(pop = F4,
                             nInd = nFamEST,
                             use = 'rand')
    F5_EST = self(F5_EST)
    F5_EST = setPheno(pop = F5_EST,
                      H2 = plotH2,
                      reps = 4, #4 Locations
                      fixEff = Fixeff,
                      p = P[b])
    # Genomic Prediction
    GS = RRBLUP(pop = F5_EST)
    
    F4_EBV = setEBV(pop = F4,
                    solution = GS,
                    value = 'bv')
    
    # Line Fixation in RGA: Season 3 (Year 2)
    F2 = self(pop = Cross,
              nProgeny = nProg)
    F3 = self(F2)
    F3 = setPheno(pop = F3,
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
    Mean$meanG_Sg1[iCycle] =  meanG(F4_EBV)
    Mean$varG_Sg1[iCycle] = varG(F4_EBV)[1]
    Mean$Accu[iCycle] = cor(gv(F4_EBV),ebv(F4_EBV))
    Mean$Cycle[iCycle] = iCycle-1
    # Stage 1 (ESTs)
    Herr$EST[b] = varG(F5_EST)[1]/varP(F5_EST)[1]
   
    #Stage 2 (AYTs)
    Mean$meanG_Sg2[iCycle] =  meanG(F6_AYT)
    Mean$varG_Sg2[iCycle] = varG(F6_AYT)[1]
    Herr$AYT[b] = varG(F6_AYT)[1]/varP(F6_AYT)[1]
   
    #AllFreq
    qtl[b+1,] = allFreq(pullQtlGeno(newParents))$maf
    alphaP[b+1,] = alphaFreq(pullQtlGeno(newParents))$p
    snp[b+1,] = allFreq(pullSnpGeno(newParents))$maf

    tmp$Scenario[b+1] = scenarioName[4]
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
          paste0("RES_",scenarioName[4], "_IRRI-OneRice_",g,"_",m,".rds"))
  
}#END OF ITERATIONS
