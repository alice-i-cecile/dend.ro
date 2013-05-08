# Core GSS functions #####

# Simulate trees under a GSS given their years of birth
sim_GSS <- function (births, last_year, GSS){
  n <- length (births)
  
  # Draw a Q value for the tree
  tree_Q <- GSS$Q0_r(n)
  
  
  # Find the natural date of death for each tree
  deaths <- births
  
  # Draw survivorship score for each tree
  score <- runif(n)
  
  # Draw lifespan of each tree
  for (i in 1:n){
    max_age <- last_year - births[i]
    Q <- tree_Q[i]
    
    # Check if the tree is surviving at the end of the simulation    
    if (score [i] <= GSS$surv_func(Q, max_age)){
      lifespan <- max_age
    } else {
      
      # Otherwise, find the year in which the fraction of the population given by the survivorship score died      
      surv_search <- function (age){
        out <- GSS$surv_func (Q, age) - score [i]
        
        return (out)
      }
      
      lifespan <- uniroot (surv_search, c(0, max_age))$root
      
      # If the tree died between years, truncate it to the last year before its death as the next ring wouldn't be fully formed      
      lifespan <- floor (lifespan)
      
    }
    
    # Year of death is year of birth + life span
    deaths[i] <- births[i] + lifespan
  }
  
  # Compile outputs into a data.frame
  demographics <- data.frame(cv_Q=tree_Q, births=births, deaths=deaths)
  
  return (demographics)
  
}

# Simulate a population under a GSS with modern sampling and uniform birth rate
modernGSS <- function (n, first_year, last_year, GSS, step_size=round(10*n)){
  
  # Set up a uniform birth distribution
  uniform_birth_dist <- function (size){
    return(sample(first_year:last_year, size, T))
  }
  
  
  # Loop structures
  enough_trees <- FALSE
  all_trees <- data.frame(cv_Q=NA, births=NA, deaths=NA)[0,]
  
  while(!enough_trees){
    
    # Decide when the trees will be born
    new_births <- uniform_birth_dist(step_size)
    
    # Draw some trees
    new_trees <- sim_GSS (new_births, last_year, GSS)
    
    # Only select the ones that are alive at the time of sampling
    surviving_trees <- new_trees[new_trees$death==last_year,]
    
    # Add them to the list
    all_trees <- rbind (all_trees, surviving_trees)
    
    # Check if there's enough
    if (nrow(all_trees)>=n){
      enough_trees <- TRUE
    }
    
  }
  
  # Once there's enough, sample at random to reach the desired number
  sampled_trees <- all_trees[sample(1:nrow(all_trees), n),]
  
  # Clean up the names
  rownames(sampled_trees) <- paste("Tree", 1:n, sep=".")

  return (sampled_trees)
}

# GSS plotting functions ####

# Year-over-year survival for plotting
getSurvProb <- function (survFunc){
  out <- function (age, growth){
    # Chance of surviving that year. (Equal to ratio of population sizes)
    stepSurvOdds <- survFunc (growth, age)/survFunc(growth, age-1)
    
    return (stepSurvOdds)
  }
  return (out)
}

# Constructing survival functions ####

make_surv_func <- function (lambda=1, k=1, Beta=0){

# lambda is shape
# k is scale
# Beta is skew for Q 

	surv_func <- function (Q, age){
	
		# Modified cumulative Weibull

	#   
	#   # Set Beta > 0 to favour slow growing trees, < 0 to favour fast growing trees
	#   # Pivots line around the point (1,1)  conserve scale at Q = 1

		Qfactor <- Beta *(Q-1) + 1

		qScale <- k*Qfactor
		
		survival <-  exp(- (age/qScale)^lambda)

		return (survival)
	}
	
	return (surv_func)

}


make_Q_0_r <- function (s){

	Q_0 <- function (x){
		return (rlnorm (x, 0, s)) 
	}
	
	return (Q_0)

}

make_Q_0_d <- function (s){

	Q_0 <- function (x){
		return (dlnorm (x, 0, s)) 
	}
	
	return (Q_0)

}


# Compile the GSS ####

make_GSS <- function (s, lambda=1, k=1, Beta=0) {
	Q0_r <- make_Q_0_r(s)
	Q0_d <- make_Q_0_d(s)
	surv_func <- make_surv_func(lambda,k,Beta)
	GSS <- list(Q0_r=Q0_r, Q0_d=Q0_d, surv_func=surv_func)
	
	return (GSS)	
}
