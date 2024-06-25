
# Allele frequencies calculation ------------------------------------------
allFreq <- function(geno) {
  n0 = apply(geno == 0, MARGIN = 2, FUN = sum) # total individuals not carrying the referent allele
  n1 = apply(geno == 1, MARGIN = 2, FUN = sum)# total individuals carrying 1 copy of the referent allele
  n2 = apply(geno == 2, MARGIN = 2, FUN = sum)# total individuals carrying 2 copies of the referent allele
  N = n0 + n1 + n2 # equal to the sample size if there are not missing data
  p = (2 * n2 + n1) / (2 * N)
  q = 1 - p
  maf = pmin(p, q) #minor allele frequency (MAF) the second most frequent allele value
  return(list(
    maf = maf,
    p = p,
    q = q,
    sizePop = unique(N)
  ))
}

# Positive allele frequencies ------------------------------------------
alphaFreq <- function(qtl) {
  addef = SP[["traits"]][[1]]@addEff
  a = NULL
  for (i in addef) {
    if (i > 0){
      ai = 2
    }else{
      ai = 0
    }
    a = append(a,ai)
  }
  
  n = NULL
  for (i in 1:nrow(qtl)) {
    tmp = qtl[i,]==a
    n = dplyr::bind_rows(n,tmp)
  }
  n = colSums(n)
  n1 = apply(qtl == 1, MARGIN = 2, FUN = sum)
  N = 2*nrow(qtl)
  p = (2*n+n1)/N
  q = 1-p
  maf = pmin(p, q)
  return(list(
    maf = maf,
    p = p,
    q = q,
    sizePop = unique(N)
  ))
}

# EBVs Prediction (rrBPLUP) -------------------------------
#require(rrBLUP)
EBV <- function(pop) {
  M = pullSnpGeno(pop) # design matrix (n × m)
  M = M - 1
  y = pop@pheno
  K = rrBLUP::A.mat(M)
  #predict marker effects
  ans = rrBLUP::mixed.solve(y, K = K)
  return(as.matrix(ans$u))
}

# EBVs Prediction (sommer) -------------------------------
#require(sommer)
EBV2 <- function(pop) {
  M = pullSnpGeno(pop) # design matrix (n × m)
  M = M - 1
  K = sommer::A.mat(M)
  dataset = matrix(nrow = pop@nInd, ncol = 2)
  colnames(dataset) = c("gid", "Trait")
  dataset = as.data.frame(dataset)
  dataset$gid = as.factor(pop@id)
  dataset$Trait = pop@pheno
  #predict marker effects
  ans <- sommer::mmer(
    fixed = Trait ~ 1,
    random = ~  vsr(gid, Gu = K),
    rcov = ~ units,
    data = dataset,
    verbose = FALSE)
  gebv= as.matrix(ans[["U"]][["u:gid"]][["Trait"]])
  return(gebv)
}
# Function for averaging multiple locations ----------------------------------
# Averaging multiple locations is done by calculating a p-value
# for the average value of the environmental covariate
pMu = function(p){
  pnorm(mean(qnorm(p)))
}

#Compute Ne from LD -------------------------------
ne.LD <- function(pop, filter = TRUE, pVal = 0.05, nChr = 12) {
  
}

#Simulation rune time -------------------------------
runTime <- function(time1, time2) {
  diffTime = as.numeric(difftime(time2, time1, units = "days"))
  #Convert the time into days/hours/mins/secs
  days = diffTime %/% 1
  Tleft = (diffTime - days) * 24
  hours = Tleft %/% 1
  Tleft = (Tleft - hours) * 60
  mins = Tleft %/% 1
  Tleft = (Tleft - mins) * 60
  secs = ceiling(Tleft)
  
  if (days >= 1) {
    Time = paste0(days, "-", hours, ":", mins, ":", secs)
  } else if (hours >= 1) {
    Time = paste0(hours, ":", mins, ":", secs)
  } else if (mins >= 1) {
    Time = paste0(mins, "m ", secs, "s")
  } else{
    Time = paste0(secs, " secs")
  }
  return(Time)
}
