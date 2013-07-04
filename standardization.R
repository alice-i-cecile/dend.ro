
# Wrapper for standardizing tree ring data ####
# tra: the tree-ring array data structure containing the data to be analysed
# effects: which effects (individual, time, age) to include in the model?
# form: are the effects added together or multiplied together
# error: is the data drawn from a normal additive PDF or a multiplicative log-normal PDF
# method: which algorithm should we use to standardize the data?

standardize_tra <- function(tra, effects=c(I=FALSE, T=TRUE, A=TRUE), form="additive", error="lnorm", method="likelihood", ...)
{
  if (method=="likelihood")
  {
    out <- standardize_likelihood(tra, effects, form, error, ...)
  }
  else if(method == "least_squares")
  {
    out <- standardize_least_squares(tra, effects, form, error, ...)
  }
  else if(method == "sfs")
  {
    out <- standardize_sfs(tra, effects, form, error, ...)
  }
  else if(method == "rcs")
  {
    out <- standardize_rcs(tra, effects, form, error, ...)
  }
  return (out)
}

# Utility functions ####

# Geometric mean function
geomMean <- function(x){
  if (length(x)==0){
    return(NA)
  }
  val <- x[!is.na(x)]
  l_val <- log(val)
  out <- exp(mean(l_val))
  return (out)
}

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
naive_effects <- function (tra, factor.dim, mean_type="arithmetic"){  
  if (mean_type=="geometric"){
    est_effect <- apply(tra, factor.dim, geomMean) 
  } else {
    est_effect <- apply(tra, factor.dim, mean, na.rm=TRUE)
  }
  
  return (est_effect)
}

# Remove effects of any number of canonical vectors from a tree ring array
remove_effect <- function (tra, effect, factor.dim, form="multiplicative"){
  
  FUN <- ifelse(form =="additive","-", "/")
  
  removed.tra <- sweep(tra, factor.dim, effect, FUN)
  
  return (removed.tra)
  
}

# Rescale effect vectors to canonical form
rescale_effects <- function (cv, form="multiplicative"){
  # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
  # So log transform, set mean to 0, then unlog
  if (form="multiplicative")
  {
    effects <- log(effects)
  }
  
  mean_effects <- lapply(effects, mean, na.rm=T)  
  
  # Scale I and T to mean of 0
  # Scale A so sum of effects stays the same
  
  # If A is missing, leave effects uncscaled
  if(!sum(!is.null(effects[[3]])))
  {
    return (effects)
  }
  
  # I
  if (sum (!is.null(effects[[1]]))){
    effects[[1]] <- effects[[1]]-mean_effects[[1]]
    effects[[3]] <- effects[[3]]+mean_effects[[1]]
  }
  
  # T
  if (sum (!is.null(effects[[2]]))){
    effects[[2]] <- effects[[2]]-mean_effects[[2]]
    effects[[3]] <- effects[[3]]+mean_effects[[2]]
  }
  
  if (form="multiplicative")
  {
    scaled_effects <- exp(scaled_effects)
  }
  
  return (effects)
}

pad_cv <- function(cv, tra, form="multiplicative"){
  
  # Set the value to fill dummy coefficients with
  if (form=="multiplicative"){
    na.value <- 1
  } else {
    na.value <- 0
  }
  
  # Initialize dummy cv lists
  new_cv <- list(Q=NA, F=NA, A=NA)
  
  # Fill empty values  
  for(i in c("Q", "F", "A")){
    if (length(cv[[i]] > 0)){
      new_cv[[i]] <- cv[[i]]
    } else {
      tra_dim <- which(c("Q", "F", "A")==i)
      
      new_cv[[i]] <- rep.int(na.value, dim(tra)[tra_dim])
      names(new_cv[[i]]) <- dimnames(tra)[[tra_dim]]
    }
  }
  
  return(new_cv)
}

sort_cv <- function(cv, tra){
  foreach(i=1:3) %do% {
    correct_order <- dimnames(tra)[[i]]
    
    cv[[i]] <- cv[[i]][correct_order]
    
  }
  
  return(cv)
}


# Regional curve standardization ####
# effect_order: the order in which effects are sequentially estimated
standardize_rcs <- function(tra, effects=c(I=FALSE, T=TRUE, A=TRUE), form="additive", error="lnorm", effect_order=c(3,2))
{
  # Select appropriate type of mean
  if (error="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Make a dummy list of effects
  effects <- list(Q=NULL, F=NULL, A=NULL)
  
  # Estimate the effects one at a time
  for (effect in effect_order){
    # Estimate an effect    
    effects[[effect]] <- naive_effect(tra, effect, mean_type)
    
    # Remove the effect
    tra <- remove_effect(tra, effects[[effect]], effect, form)
  }
  
  # Fill in dummy values for effects not estimated
  cv <- pad_cv(cv, tra, form)
  
  # Make sure effects are in the right order
  cv <- sort_cv(cv, tra)
  
  # Rescale the CV to standard form
  cv <- rescaleCV(cv, form)
  
  names (cv) <- c("I","T","A")
  
  return (list(cv=cv))   
}