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

# Get sample size along a dimension ####
sample_depth_tra <- function(tra, factor.dim=2){ #1 is tree, 2 is time, 3 is age
  
  filled_cells <- !is.na(tra)
  
  sample_depth <- apply(filled_cells, factor.dim, sum)

  return (sample_depth)
}
