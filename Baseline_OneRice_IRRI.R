m = Sys.getenv("SLURM_ARRAY_TASK_ID")
m = as.integer(m)

nReps=1

for (REP in 1:nReps) {
  
  load(paste0("BURNIN_OneRice_",g,"_",m,".RData"))
  
  cat("----------------------------------------------------\n")
  cat(paste0("RUNING OF THE BREEDING SCHEME ", scenarioName[1]," ... \n"))
  cat("----------------------------------------------------\n")
  
  for (j in 1:Cohort5Y) { 
    
    #Select parents from Base population
    Parents = selectInd(pop = ECP,
                        nInd = nParents,
                        use = "rand")
    
    # Crossing Block & F1s Validation: Season 1 & 2 (Year 1)
    Cross = randCross(pop = Parents,
                      nCrosses = nCross,
                      nProgeny = 1)
    
    if (j < 5){ 
      # Line Fixation in RGA: Season 3 & 4 (Year 2)
      F2 = self(pop = Cross,
                nProgeny = nProg)
      F3 = self(F2)
      F4 = self(F3)
      F5 = self(F4)
      F6 = self(F5)
    }
    
    if (j < 4){ 
      # Line Testing: Season 5 (Year 3)
      F6_LST = setPheno(pop = F6,
                        H2 = plantH2,
                        reps = 1, #1 Location
                        fixEff = Fixeff, 
                        p = pVal[j+2])
    }
    
    if (j < 3){ 
      # 1st Stage of Yield Evaluation (Seed Amplification): Season 7 (Year 4)
      F7_OYT = selectWithinFam(pop = F6_LST,
                               nInd = nFamOYT,
                               use = 'pheno')
      # Advancement
      F7_OYT = self(F7_OYT)
      F7_OYT = setPheno(pop = F7_OYT,
                        H2 = plantH2,
                        reps = 1, #1 Location
                        fixEff = Fixeff, 
                        p = pVal[j+3])
    }
    
    if (j < 2){ 
      # 1st Stage of Yield Evaluation (Estimation Set): Season 9 & 10 (Year 5)
      F8_EST = selectWithinFam(pop = F7_OYT,
                               nInd = nFamEST,
                               use = 'pheno')
      F8_EST = self(F8_EST)
      #F8_EST = self(F7_OYT)
      F8_EST = setPheno(pop = F8_EST,
                        H2 = plotH2, 
                        reps = 4, #4 Locations
                        fixEff = Fixeff, 
                        p = pVal[j+4])
    }
  }#END Pipeline
  
  Mean$Scenario = scenarioName[1]; Mean$Kind = "PS"
  
  # RUN 30 Years of Target Breeding Scheme
  
  for (b in iBreed:nYears) {
    
    iCycle = b-nBrn+1
    
    #Select new parents from last estimated breeding values
    newParents = selectInd(pop = F8_EST,
                           nInd = nParents,
                           use = 'pheno')
    
    # 2nd Stage of Yield Evaluation on AYT: Season 11 & 12 (Year 6)
    F9_AYT = selectInd(pop = F8_EST,
                       nInd = nAYT,
                       use = 'pheno')
    # Advancement
    F9_AYT = self(F9_AYT)
    F9_AYT = setPheno(pop = F9_AYT,
                      H2 = plotH2,
                      reps = 4, #4 Locations
                      fixEff = Fixeff, 
                      p = P[b])
    
    # 1st Stage of Yield Evaluation (Estimation Set): Season 9 & 10 (Year 5)
    F8_EST = selectWithinFam(pop = F7_OYT,
                             nInd = nFamEST,
                             use = 'pheno')
    F8_EST = self(F8_EST)
    #F8_EST = self(F7_OYT)
    F8_EST = setPheno(pop = F8_EST,
                      H2 = plotH2,
                      reps = 4, #4 Locations
                      fixEff = Fixeff,
                      p = P[b])
    
    # 1st Stage of Yield Evaluation (Seed Amplification): Season 7 (Year 4)
    F7_OYT = selectWithinFam(pop = F6_LST,
                             nInd = nFamOYT,
                             use = 'pheno')
    # Advancement
    F7_OYT = self(F7_OYT)
    F7_OYT = setPheno(pop = F7_OYT,
                      H2 = plantH2,
                      reps = 1, #1 Location
                      fixEff = Fixeff, 
                      p = P[b])
    
    # Line Testing: Season 5 (Year 3)
    F6_LST = setPheno(pop = F6,
                      H2 = plantH2,
                      reps = 1, #1 Location
                      fixEff = Fixeff, 
                      p = P[b])
    
    # Line Fixation in RGA: Season 3 & 4 (Year 2)
    F2 = self(pop = Cross,
              nProgeny = nProg)
    F3 = self(F2)
    F4 = self(F3)
    F5 = self(F4)
    F6 = self(F5)
    
    # Crossing Block & F1s Validation: Season 1 & 2 (Year 1)
    Cross = randCross(pop = newParents,
                      nCrosses = nCross,
                      nProgeny = 1)
    
    # Save results
    cat(paste0("Saving Data From Crossing Year ", b,  " ....\n"))
    
    #LSTs
    Herr$LST[b] = varG(F6_LST)[1]/varP(F6_LST)[1]
    #OYTs
    Herr$OYT[b] = varG(F7_OYT)[1]/varP(F7_OYT)[1]
    
    # Stage 1 (ESTs)
    Mean$meanG_Sg1[iCycle] =  meanG(F8_EST)
    Mean$varG_Sg1[iCycle] = varG(F8_EST)[1]
    Mean$Accu[iCycle] = cor(gv(F8_EST),EBV(F8_EST))
    Mean$Cycle[iCycle] = iCycle-1
    Herr$EST[b] = varG(F8_EST)[1]/varP(F8_EST)[1]
    
    #Stage 2 (AYTs)
    Mean$meanG_Sg2[iCycle] =  meanG(F9_AYT)
    Mean$varG_Sg2[iCycle] = varG(F9_AYT)[1]
    Herr$AYT[b] = varG(F9_AYT)[1]/varP(F9_AYT)[1]
    
    #AllFreq
    qtl[b+1,] = allFreq(pullQtlGeno(newParents))$maf
    alphaP[b+1,] = alphaFreq(pullQtlGeno(newParents))$p
    snp[b+1,] = allFreq(pullSnpGeno(newParents))$maf
    
    tmp$Scenario[b+1] = scenarioName[1]
    tmp$Kind[b+1] = "PS"
    
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
          paste0("RES_",scenarioName[1], "_IRRI-OneRice_",g,"_",m,".rds"))
  
}#END OF ITERATIONS

warnings()