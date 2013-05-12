# Core functions ####

# Naive estimate of canonical vectors (analagous to constructing regional curve)
naiveCV <- function (tra, factor.dim, meanType="arithmetic"){
# Traditional RCS uses the arithmetic mean, but should use the geometric mean
  
  if (meanType=="geometric"){
    est_cv <- apply(tra, factor.dim, geomMean) 
  } else {
    est_cv <- apply(tra, factor.dim, mean, na.rm=TRUE)
  }
  
  return (est_cv)
}

# Remove effects of any number of canonical vectors from a tree ring array
removeCV <- function (tra, cv, factor.dim, FUN="/"){
  
  removed.tra <- sweep(tra, factor.dim, cv, FUN)
  
  return (removed.tra)
  
}

# Regional curve standardization ####
standardize_rcs <- function (tra, factor_order=c(3,2), meanType="arithmetic") {
  if (meanType=="arithmetic"){
    multiplicative <- FALSE
  } else {
    multiplicative <- TRUE
  }
  
  cv <- list(Q=NULL, F=NULL, A=NULL)
  

  # Estimate the effects one at a time
  for (effect in factor_order){
    # Estimate an effect    
    cv[[effect]] <- naiveCV(tra, effect, meanType)
    
    # Remove the effect
    tra <- removeCV(tra, cv[[effect]], effect)
    
  }
  
  # Fill in dummy values for effects not estimated
  cv <- pad_cv(cv, tra, multiplicative=FALSE)
  
  # Make sure effects are in the right order
  cv <- sort_cv(cv, tra)
  
  # Rescale the CV to standard form
  cv <- rescaleCV(cv)
  
  names (cv) <- c("Q","F","A")
  
  return (list(cv=cv))  
}

# Signal-free regional curve standardization ####
standardize_sfs <- function (tra, factor_order=c(3,2), cor_threshold=0.999999, meanType="arithmetic"){
  
  # Factor order needs to be of length 2 to work  
  
  converged <- FALSE
  iteration <- 0
  
  working_tra <- tra
  
  while (!converged){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
        
    # Estimate the values for the first dimension
    est_1 <- naiveCV(working_tra, factor_order[1], meanType)
    #print(est_1)
          
    # Remove those effects temporarily
    intermediate_tra <- removeCV (working_tra, est_1, factor_order[1])
    
    # Estimate values for the second dimensions
    est_2 <- naiveCV(intermediate_tra, factor_order[2], meanType)
    #print(est_2)
     
    # Remove those effects from the working data
    working_tra <- removeCV (working_tra, est_2, factor_order[2])
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal
    if (meanType=="arithmetic"){
      conv_cor <- cor(working_tra, last_tra,  "complete.obs") 
    } else {
      conv_cor <- cor(log(working_tra), log(last_tra),  "complete.obs")
    }
    
    print (paste("Log-correlation of current and last iteration of", conv_cor))
  
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
    
  }
  
  # Create storage for the estimated CV
  all_cv <-vector(mode="list", length=3)

  # Dummy missing CV
  all.dim <- c(1,2,3)
  
  miss.dim <- all.dim[-intersect(factor_order, all.dim)]
  miss.cv <- rep.int(1, dim(tra)[[miss.dim]])
  names (miss.cv) <- dimnames(tra)[[miss.dim]]
  
  # Primary chronology is mean of converged working TRA
  prim.cv <- naiveCV(working_tra, factor_order[1], meanType)
  
  # Secondary chronology is mean of original data with primary effects removed
  sec.series <- removeCV (tra, prim.cv, factor_order[1])
  sec.cv <- naiveCV(sec.series, factor_order[2], meanType)
  
  # Compile the CV in the approriate order
  all_cv <- list(prim.cv, sec.cv, miss.cv)
  
  cv.order <- order(c(factor_order[1], factor_order[2], miss.dim))
  all_cv <- all_cv[cv.order]
    
  # Rescale the CV to standard form
  all_cv <- rescaleCV(all_cv)
  
  return (list(cv=all_cv))
}

# Truly signal-free regional curve standardization ####
# Cleans up SF-RCS algorithm and allows expansion to N dimensions

standardize_tsfs <- function (tra, factor_order=c(3,2), cor_threshold=0.999999, meanType="arithmetic", corr_m=T){
  
  # Total dimensionality of the data
  tra_dim <- length(dim(tra))
  
  # Create storage for the estimated CV
  all_cv <-vector(mode="list", length=tra_dim)
  
  # Dummy starting CVs
  for (i in 1:tra_dim){
    all_cv[[i]] <-  rep.int(1, dim(tra)[[i]])
    names(all_cv[[i]]) <- dimnames(tra)[[i]]
  }
  
  # Loop controls
  converged <- FALSE
  iteration <- 0
  working_tra <- tra
  
  while (!converged){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
    
    for (j in factor_order){
      
      # Estimate the effects across each dimension
      est_j <- naiveCV(working_tra, j, meanType)
     
      # Remove them from the signal-free data
      working_tra <- removeCV (working_tra, est_j, j)
      
      # Combine them with previously determined effects for that dimension
      all_cv[[j]] <- all_cv[[j]] * est_j
      #print(all_cv[[j]])
    }
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal
    if (meanType=="arithmetic"){
      conv_cor <- cor(working_tra, last_tra,  "complete.obs") 
    } else {
      conv_cor <- cor(log(working_tra), log(last_tra),  "complete.obs")
    }
    
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
    
  }
  
  # Rescale the CV to standard form
  all_cv <- rescaleCV(all_cv)
  
  names (all_cv) <- c("Q","F","A")

  return (list(cv=all_cv))
}
