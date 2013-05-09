# Dependencies ####
source ("GSS.R")

# Trends for age ####

# Trend for constant basal area increment
constBAI_trend <- function (x, k, p0=0){
  # c is the constant basal area increment
  # p0 is the pith offset
  out <- sqrt(k*(x)/pi+p0^2) - sqrt(k*(x-1)/pi+p0^2)
  return (out)
}

# Trend for flat growth with time
flat_trend <- function(x, k=1){
  return (k)
}

linear_trend <- function(x, k=1, s=0){
  return (k+s*x)
}

salty_trend <- function(x, k=1){
  return(k*rlnorm(1,0,0.01))
}

# Trend for negative exponential growth with time
exp_trend <- function(x, k=-1, t=0){
  return (exp(-k*x)+t)
}

rand_trend <- function (x, noise=0.1){
  return (rlnorm(1, 0, noise))
}

rand_lnorm_trend <- function (x, m=1, s=0){
  return (rlnorm(1, m, s))
}

# Adding noise to synthetic TRAs ####
add_noise <- function (tra, noiseSD, multiplicative=TRUE){
  noisy_tra <- tra
  
  if (multiplicative){
    noisy_tra [!is.na(tra)] <- tra [!is.na(tra)]*rlnorm(sum(!is.na(tra)), 0, noiseSD)   
  } else {
    noisy_tra [!is.na(tra)] <- tra [!is.na(tra)]+rnorm(sum(!is.na(tra)), 0, noiseSD)
    
    # Negative and zero values are invalid, remove them
    noisy_tra[noisy_tra <= 0] <- NA
  }
  
  return (noisy_tra)
}

# Generate simple synthetic TRAs ####
crude_synth_TRA <- function (nQ, nF, nA, funcA, meanQ=0, meanF=0, meanA=1, sdQ=0, sdF=0, sdA=0, noiseSD=0, retentionFraction=1, multiplicative=T, ...){
  
  # Random CVs
  cv.Q <- rlnorm (nQ, meanQ, sdQ)
  names(cv.Q) <- paste ("Q", 1:nQ, sep="")
  cv.F <- rlnorm (nF, meanF, sdF)
  names(cv.F) <-  paste ("F", 1:nF, sep="")
  #cv.A <- rlnorm (nA, meanA, sdA)
  #cv.A <- sapply(1:nA, funcA, m=meanA, s=sdA)
  cv.A <- sapply(1:nA, funcA, ...)
  names(cv.A) <-  paste ("A", 1:nA, sep="")
  
  cv <- list("Q"=cv.Q, "F"=cv.F, "A"=cv.A)
  
  cv <- rescaleCV(cv)
  
  tra <- cv[[1]] %o% cv[[2]] %o% cv[[3]]
  
  # Blank out a fraction of the data at random
  weight_list <- rbinom(length(tra), 1, retentionFraction)
  weight_list [weight_list==0] <- NA
  
  tra <- tra * weight_list
  
  # Add  noise
  tra <- add_noise (tra, noiseSD, multiplicative)
  
  return (list(cv=cv, tra=tra))
}

# Realistic TRAs ####

# Base constructor function
base_synth_TRA <- function (cv_Q, cv_F, cv_A, births, deaths, noiseSD=0, multiplicative=TRUE){
  
  # Make a generic simple RWL for the appropriate birth and death years
  rwl <- as.data.frame(matrix(1, length(cv_F), length(cv_Q), dimnames=list(names(cv_F), names(cv_Q))))
  rownames(rwl) <- names(cv_F)
    
  c_years <- as.numeric(rownames(rwl))  
  for (i in 1:length(rwl)){
    # Find years the tree was alive
    living_years <- c_years >= births[i] & c_years <= deaths[i]
        
    # Set all values outside of living years to NA
    # Living years have a base value of 1
    rwl[!living_years, i] <- NA    
    
  }
    
  # Convert to a TRA
  tra <- rwl.to.tra(rwl)
    
  # Truncate the CV according to the values actually observed
  cv_F <- cv_F[dimnames(tra)[[2]]]
  cv_A <- cv_A[dimnames(tra)[[3]]]
    
  # Scale the TRA according to the CV
  tra <- sweep(tra, 1, cv_Q, "*")
  tra <- sweep(tra, 2, cv_F, "*")
  tra <- sweep(tra, 3, cv_A, "*")
  
  # Add noise
  tra <- add_noise (tra, noiseSD, multiplicative)
  
  return (tra)
}

# Basic realistic parameterization and construction for a modern data set

modern_TRA <- function (nQ, nF, sdF, funcA, noiseSD=0, multiplicative=TRUE, GSS, ...){
  # Assumes:
  # constant birth rate
  # Log-normal random forcing and noise
  # Strictly modern sampling
  
  # Generate Q, birth and death dates using the GSS
  demographics <- modernGSS(nQ, 1, nF, GSS)
  
  cv_Q <- demographics$cv_Q
  #names(cv_Q) <- rownames(demographics)
  names(cv_Q) < 1:nQ
  
  births <- demographics$births
  deaths <- demographics$deaths
  
  # Generate forcing
  cv_F <- rlnorm(nF, 0, sdF)
  names(cv_F) <- 1:nF
  
  
  # Generate age trend
  cv_A <- sapply(1:nF, funcA, ...)
  names(cv_A) <- 1:nF
    
  # Construct full TRA
  tra <- base_synth_TRA(cv_Q, cv_F, cv_A, births, deaths, noiseSD, multiplicative)  
  
  # Truncate age trend according to oldest tree observed
  cv_A <- cv_A[1:dim(tra)[[3]]]
  
  # Bundle CVs  
  cv <- list("Q"=cv_Q, "F"=cv_F, "A"=cv_A)
  
  # Rescale CVs
  cv <- rescaleCV(cv)
  
  # Compile output
  out <- list(tra=tra, cv=cv)
  
  return(out)
}

correlate_TRA <- function (tra, cv, corr){
  
  ageList <- rep.int(NA, length(cv[[1]]))
  names(ageList) <- names (cv[[1]])
  
  for (i in 1:length(cv[[1]])){
    ageList[i] <- max(which (!is.na(tra[i,,]), arr.ind=T)[,2])
  }

  # Generate new correlated Q values
#   new_Q <- corr*ageList + sqrt(1-corr^2)*log(cv[[1]])
    new_Q <-corr*ageList + sqrt(1-corr^2)*cv[[1]]

  # Rescale them so geomMean = 1
#   new_Q <- new_Q-mean(new_Q)
    new_Q <- new_Q/geomMean(new_Q)
  
  # Untransform
#   new_Q <- exp(new_Q)
  
  # Recombine
  new_cv <- list(new_Q, cv[[2]], cv[[3]])
  
  # Construct new data
  new_tra <- new_cv[[1]] %o% new_cv[[2]] %o% new_cv[[3]]
  
  # Cull to match form
  new_tra [is.na(tra)] <- NA
  
  
  return (list(tra=new_tra, cv=new_cv))
  
}
