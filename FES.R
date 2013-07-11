# Load dependencies ####
library(foreach)
library(gam)

# Core regression function ####

# Perform factor regression standardization on a tree ring array
standardize_fes <- function (tra, incQ=T, incF=T, incA=T, multiplicative=T, corr_m=T, model_type="lm", span="auto", degree=1, ...){
  
  
  # Sanitize 0 values to NA since ringwidth of zero is really a missing value
  tra [which(tra==0)] <- NA
  
  if (multiplicative){
    # Log-transform the data
    tra <- log(tra)
  }
  
  # Convert into a useable format
  tab.tra <- tra.to.table (tra)
  
  # Sanitize data frames so Inf / NaN values are recorded as NA
  tab.tra[!is.finite(tab.tra$G),"G"] <- NA
  
  # Construct formula for regression
  if (model_type=="lm"){
    growth_formula <- as.formula(make_lm_formula(incQ, incF, incA))
  } else if (model_type=="gnm") {  
    growth_formula <- as.formula(make_gnm_formula(incQ, incF, incA))
  } else if (model_type=="gam"){
    growth_formula <- as.formula(make_gam_formula(incQ, incF, incA))
  }
  
  print ("Model initialized.")
  print (growth_formula)
  
  # Find the canonical vectors using regression
  
  # If error is multiplicative log-normal:  
  # Use linear regression on the log-transformed data
  # lm, biglm, rlm, glm also work with same sytax
  # gnm can be used to provide common interface for the two approaches
  
  # If error is additive and normal:
  # Use nonlinear regression on original model and data
  # Much slower, gnm used for ease of working with factors
  
  if (model_type=="lm"){
    # Linear fit
    growth_model <- lm (growth_formula, data=tab.tra,...)
  } else if (model_type=="gnm"){
    # Nonlinear fit
    growth_model <- gnm (growth_formula, data=tab.tra, ...)
  } else if (model_type=="gam"){
    
    # GAM FES is equivalent to the faster LM FES if the age effect is not computed
    if(incA==FALSE){
      growth_model <- lm (growth_formula, data=tab.tra,...)
    } else {
      # Make a numeric version of the age factor for smoothing
      tab.tra$A.num <- as.numeric(as.character(tab.tra$A.cv))
      
      # Minimum span is (degree+1)/Number of distinct observations
      # Maximum span is 1
      min_span <- (degree+1)/length(unique(tab.tra$A.num))
      max_span <- 1
      
      # If no smoothing span parameter is given, select one automatically using AIC
      if (span=="auto"){
        
        print ("Automatically selecting span parameter by AIC.")
        
        optim_span <- function (span){
          local_formula <- growth_formula
          environment (local_formula) <- environment()
          model <- gam(local_formula, data=tab.tra)
          
          model_AIC <- AIC(model)
          print(c(Span=span, AIC=model_AIC))        
          return(model_AIC)
        }
        
        span <- optimize(f=optim_span, interval=c(min_span, max_span))$minimum
        
        print (paste("Span of", span, "selected."))
        
      } else {
        if (span < min_span){
          span <- min_span
          print(paste("Selected span too low, automatically set to", min_span))
        }
        if (span > max_span){
          span <- max_span
          print(paste("Selected span too high, automatically set to", 1))
        }
      } 
      
      # GAM smoothed fit
      growth_model <- gam (growth_formula, data=tab.tra, ...)
    }
  }
  
  model_summ <- summary(growth_model)
  print ("Model fit complete.")
  
  
  # Extract the coefficients found
  cv <- extract_cv (growth_model, model_type)
  cv <- pad_cv (cv, tra, multiplicative)  
  cv <- sort_cv (cv, tra)
        
  print("Canonical vectors and standard errors of estimate extracted.")
  
  if (multiplicative){
    # Untransform the canonical vectors
    full_cv <- lapply(cv, length)>0    
    cv[full_cv] <- lapply(cv[full_cv], exp) 
  }
  
  # Scale the canonical vectors and confidence intervals
  cv <- rescaleCV (cv)
  
  # Compute the standard deviation of the residuals
  sd.epsilon <- sd(residuals(growth_model), na.rm=T)
  
  print (paste("Standard deviation of residuals:", signif(sd.epsilon,3))) 

  # Compute the R^2 values for the fit  
  if (model_type=="lm"){
    Rsq <-  model_summ$r.squared
    adj.Rsq <- model_summ$adj.r.squared
  }

  # Compute the log-likelihood of the model
  loglike <- logLik (growth_model)
  
  # Compute the AIC of the fit
  aic <- AIC(growth_model)
  bic <- BIC(growth_model)
  
  # Compile outputs  
  if (model_type=="lm"){
    fit <- c(Rsq, adj.Rsq, loglike, aic, bic)
    names (fit) <- c("R-squared", "Adjusted.R-squared", "Log-Likelihood", "AIC", "BIC")
  } else {
    fit <- c(loglike, aic, bic)
    names (fit) <- c("Log-Likelihood", "AIC", "BIC")
  }
  
  print ("Fit computed.")
  print (fit)
  
  names (cv) <- c("Q","F","A")
    
  # Standard errors returned are typically log-transformed, should be converted before plotting
  output <- list("cv"=cv, "RSE"=sd.epsilon, "fit"=fit, "model"=growth_model)
  
  return (output)
}

# Tree ring tables ####

# Convert a tree ring array to a table suitable for use in lm () type functions
tra.to.table <- function (tra){
  
  full <- !is.na(tra)
  values <- tra[full]
  ind_pos <- which(full, arr.ind=T)
  
  
  Q_pos <- dimnames(full)[[1]][ind_pos[,1]]
  F_pos <- dimnames(full)[[2]][ind_pos[,2]]
  A_pos <- dimnames(full)[[3]][ind_pos[,3]]
  
  df <- data.frame(G=values, Q.cv=Q_pos, F.cv=F_pos, A.cv=A_pos)
  
  return (df)
}

# Formula construction ####
make_gam_formula <- function (incQ, incF, incA){
  dep.str <- "G"
  ind.str <- "0"
  
  if(incQ){
    ind.str <- paste(ind.str, "Q.cv", sep="+")
  }
  if(incF){
    ind.str <- paste(ind.str, "F.cv", sep="+")
  }
  if(incA){
    # Change smoothing terms here
    ind.str <- paste(ind.str, "lo(A.num, span=span, degree=degree)", sep="+")
  }
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}

make_lm_formula <- function (incQ, incF, incA){
  dep.str <- "G"
  ind.str <- "0"
  
  if(incQ){
    ind.str <- paste(ind.str, "Q.cv", sep="+")
  }
  if(incF){
    ind.str <- paste(ind.str, "F.cv", sep="+")
  }
  if(incA){
    ind.str <- paste(ind.str, "A.cv", sep="+")
  }
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}

make_gnm_formula  <- function(incQ, incF, incA){
  # Determine formula desired
  dep.str <- "G"
  
  # Combine the independent variables    
  vars.str <- ""
  
  if(incQ){
    vars.str <- paste(vars.str, "Q.cv", sep=",")
  }
  if(incF){
    vars.str <- paste(vars.str, "F.cv", sep=",")
  }
  if(incA){
    vars.str <- paste(vars.str, "A.cv", sep=",")
  }
  
  # Remove leading symbol
  vars.str <- substr(vars.str, 2, nchar(vars.str))
  
  # Add in the Mult() term and set the offset to 0
  ind.str <- paste("Mult(", vars.str, ")+0", sep="")
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}

# Extracting info from a model fit ####
extract_cv <- function (growth_model, model_type){
  model_class <- class(growth_model)
  
  if (model_type=="gam"){
    cv <- extract_cv_gam(growth_model)
  } else if (model_type=="gnm"){
    cv <- extract_cv_gnm(growth_model)
  } else if (model_type=="lm"){
    cv <- extract_cv_lm(growth_model)
  }
  
  return (cv)
  
}

extract_cv_lm <- function (growth_model){
  coeff <- dummy.coef(growth_model)
  Q.cv <- coeff$Q.cv
  F.cv <- coeff$F.cv
  A.cv <- coeff$A.cv
  
  cv <- list("Q"=Q.cv, "F"=F.cv, "A"=A.cv)
  
  return (cv)
}

extract_cv_gnm <- function (growth_model){
  coeff <- coef(growth_model)
  
  names(coeff) <- gsub("Mult.*)\\.", "", names(coeff))
  Q.cv <- coeff[grep("Q.cv", names(coeff))]
  F.cv <- coeff[grep("F.cv", names(coeff))]
  A.cv <- coeff[grep("A.cv", names(coeff))]
  
  names(Q.cv) <- substr(names(Q.cv), 5, nchar(names(Q.cv)))
  names(F.cv) <- substr(names(F.cv), 5, nchar(names(F.cv)))
  names(A.cv) <- substr(names(A.cv), 5, nchar(names(A.cv)))
  
  cv <- list("Q"=Q.cv, "F"=F.cv, "A"=A.cv)
  
  # gnm function returns the negative coefficients in ~50% of the cases as they are also a least squares solution
  # Force the positive solution
  
  cv <- lapply(cv, abs)
  
  return (cv)
  
}

extract_cv_gam <- function (growth_model){  
  # Extract the linear components, Q and F
  coeff <- dummy.coef(growth_model)
  Q.cv <- coeff$Q.cv
  F.cv <- coeff$F.cv
  
  # Cross reference with original data table
  tab.tra <- growth_model$data
  smooth_fit <- growth_model$smooth
  
  unique_levels <- unique(tab.tra$A.num)
  indices <- unlist(foreach(i=unique_levels) %do% {match(i, tab.tra$A.num)})
  
  A.cv <- smooth_fit[indices]
  names(A.cv) <- unique_levels
  
  cv <- list("Q"=Q.cv, "F"=F.cv, "A"=A.cv)
  
  return (cv)
}

extract_se <- function (growth_model){
  predicted_se <- predict(growth_model, se.fit=T, type="terms")[[1]]
  
  tab.tra <- growth_model$model
  
  se <- lapply(1:ncol(predicted_se), cross_ref_model_terms, predicted_terms=predicted_se, model=tab.tra)
  
  return(se)
}
