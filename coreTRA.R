# Geometric mean utility function
geomMean <- function(x){
  if (length(x)==0){
    return(NA)
  }
  val <- x[!is.na(x)]
  l_val <- log(val)
  out <- exp(mean(l_val))
  return (out)
}

# Get endpoints of a rwl series
getEndpoints <- function (series, side="start"){
  #side: "start" or "stop" of series
  rings <- subset(series, !is.na(series))
  if (side=="start"){
    return (head(rownames(rings),1))
  }
  if (side=="stop"){
    return (tail(rownames(rings),1))
  }      
}

# Convert a standard tree ring data frame (.rwl) into tree ring array form. Model is G(i, t, T) = Q*F*A
rwl.to.tra <- function (rwl, birth_years=NULL){
  
  if (is.null(birth_years)){
    # Determine birth for each tree
    birth_years <- foreach(i=colnames(rwl)) %do% {getEndpoints(rwl[i], side="start")}
    names (birth_years) <- names(rwl)
    
  }  else {
    # If rwl file is too short, add empty rows to it
    birth_limit <- min(sapply(birth_years, as.numeric))
    rwl_limit <- min(sapply(rownames(rwl), as.numeric))
    
    if (birth_limit < rwl_limit){
      num_missing <- rwl_limit-birth_limit
      
      empty_rows <- rbind(rwl[0,], matrix(NA, num_missing, ncol(rwl)))
      
      colnames(empty_rows) <- colnames(rwl)
      rownames(empty_rows) <- birth_limit:(rwl_limit-1)
      
      rwl <- rbind(empty_rows, rwl)
    }
  }
  
  # Find the indices for each year
  birth_index <- vector ()
  
  for (i in 1:length(birth_years)){
    birth_index[i] <- which(rownames(rwl)==toString(birth_years[i]))
  }
  
  names (birth_index) <- names (birth_years)
    
  # Compute dimension size for all elements
  i.size <- ncol (rwl)
  t.size <- nrow (rwl)
  T.size <- t.size
  
#   # Construct empty vectors of the appropriate length
#   Q.empty <- rep.int (NA, i.size)
#   F.empty <- rep.int (NA, t.size)
#   A.empty <- rep.int (NA, T.size)
#   
#   # Construct an empty tree ring array via tensor multiplication
#   tra <- Q.empty %o% F.empty %o% A.empty
  
  # Construct an array directly
  tra <- array(NA, c(i.size, t.size, T.size), list(colnames(rwl), rownames (rwl), 1:T.size))
  
  # Input data 
  for (tree in 1:ncol(rwl)){
    birth <- birth_index[tree]
    for (year in 1:nrow(rwl)){
      
      # Find age of ring data
      age <- year - birth + 1
      

      datum <- rwl [year, tree]
      
      # Only update filled values
      if (!is.na(datum)){        
        tra [tree,year,age] <- datum
      }
    }
  }
  
  # Truncate empty rows
  empty_tra <- is.na(tra)
  
  empty_Q <- apply(empty_tra, 1, all)
  empty_F <- apply(empty_tra, 2, all)
  empty_A <- apply(empty_tra, 3, all)
  
  tra <- tra[!empty_Q, !empty_F, !empty_A]
  
  return (tra)
}

# Convert a tree ring array back to a standard tree ring data frame
tra.to.rwl <- function (tra) {
  
  
  # Get dimension sizes
  i.size <- dim (tra)[1]
  t.size <- dim (tra)[2]
  
  # Make an empty data frame
  rwl <- as.data.frame(matrix (NA, t.size, i.size))
  
  colnames(rwl) <- dimnames(tra)[[1]]
  rownames(rwl) <- dimnames(tra)[[2]]
  
  # Fill data
  for (tree in 1:dim(tra)[1]){
    
    # Get lifespan to fill. Could be more general in style
    filled <- !is.na(tra[tree, ,])
    filledYears <- rownames(filled[rowSums(filled)>0,])
    names(filledYears) <- dimnames(tra)[[1]]
    
    # Extract data
    data <- tra[tree,,]
    data <- data[!is.na(data)]
    
    # Put data into appropriate position in data frame
    rwl[filledYears, tree] <- data
  }
  
  return (rwl)
}

# Construct a partially filled tra from canonical vectors and a reference tree ring array describing position data
tra.from.cv <- function (cv, tra=NULL){
  
  # Construct full tree ring array
  completeTRA <-  cv[[1]]
  
  for (i in 2:length(cv)){
    completeTRA <- completeTRA %o% cv[[i]]
  }
  
  if (!is.null(tra)){
    # The position array only retains presence/absence data
    positionArray <- tra
    positionArray[!is.na(tra)] <- 1
    
    # Retain only values that were present before. Uses dimnames from the canonical vectors
    newTRA <- completeTRA * positionArray
    
    return (newTRA)
    
  } else {
    
    return (completeTRA)
  }
  
}

# Rescale canonical vectors to canonical form
rescaleCV <- function (cv, ci.u=NULL, ci.l=NULL){
  meanCV <- lapply(cv, geomMean)
  # meanCV <- lapply(cv, mean)
  
  
  # Scale Q and F to mean of 1
  # Scale A so product of canonical vectors stays the same
  
  # Q
  if (sum (!is.null(cv[[1]]))){
    cv[[1]] <- cv[[1]]/meanCV[[1]]
    cv[[3]] <- cv[[3]]*meanCV[[1]]
    
    # Also rescale confidence intervals
    if (!is.null(ci.u)){
      ci.u[[1]] <-  ci.u[[1]]/meanCV[[1]]      
      ci.u[[3]] <-  ci.u[[3]]*meanCV[[1]]
    }
    
    if (!is.null(ci.l)){
      ci.l[[1]] <-  ci.l[[1]]/meanCV[[1]]
      ci.l[[3]] <-  ci.l[[3]]*meanCV[[1]]
    }
    
  }
  
  # F
  if (sum (!is.null(cv[[2]]))>0){
    cv[[2]] <- cv[[2]]/meanCV[[2]]
    cv[[3]] <- cv[[3]]*meanCV[[2]]
    
    # Also rescale confidence intervals
    if (!is.null(ci.u)){
      ci.u[[2]] <-  ci.u[[2]]/meanCV[[2]]      
      ci.u[[3]] <-  ci.u[[3]]*meanCV[[2]]
    }
    
    if (!is.null(ci.l)){
      ci.l[[2]] <-  ci.l[[2]]/meanCV[[2]]
      ci.l[[3]] <-  ci.l[[3]]*meanCV[[2]]
    }
    
  }  
  
  if (is.null(ci.u) & is.null(ci.l)){
    return (cv)
  } else {
    
    output <- list(CV=cv, CI.U=ci.u, CI.L=ci.l)
    return (output)
  }
}

pad_cv <- function(cv, tra, multiplicative){
  
  if (multiplicative){
    na.value <- 0
  } else {
    na.value <- 1
  }
  
  new_cv <- list(NULL, NULL, NULL)
  
  compare_names <- function (names_cv, names_tra, dim){
    overlap_Q <- length(intersect(names_cv, names_tra[[1]]))
    overlap_F <- length(intersect(names_cv, names_tra[[2]]))
    overlap_A <- length(intersect(names_cv, names_tra[[3]]))
    
    if (overlap_Q==length(names_cv)){return (1)}
    if (overlap_Q==length(names_cv)){return (2)}
    if (overlap_Q==length(names_cv)){return (3)}
    
    if ((overlap_Q != overlap_F) & (overlap_Q!= overlap_A) & (overlap_F!= overlap_A)){
      best <- which.max(c(overlap_Q, overlap_F, overlap_A))
    } else {
      best <- dim
    }
        
    return (best)
  }    
  
  # Match vectors to appropriate dimensions
  for(i in 1:length(cv)){
    new_cv[[compare_names(names(cv[[i]]), dimnames(tra), i)]] <- cv[[i]]
  }
  
  # Pad empty CV with zeros / null values
  for (j in 1:3){
    if (length(new_cv[[j]])==0){
      new_cv[[j]] <- rep.int(na.value, dim(tra)[j])
      names(new_cv[[j]]) <- dimnames(tra)[[j]]
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



# Get sample size along a dimension ####
sample_depth_tra <- function(tra, factor.dim=2){ #1 is tree, 2 is time, 3 is age
  
  filled_cells <- !is.na(tra)
  
  sample_depth <- apply(filled_cells, factor.dim, sum)

  return (sample_depth)
}
