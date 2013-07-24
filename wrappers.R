# Load dependencies ####
library (stringr)
library (dplR)  
library (MASS)
library (reshape)
library (car)
library (nlme)
library (gnm)
library (gtools)
library (plyr)
library (rjson)

source ("coreTRA.R")
source ("standardization.R")
source ("GSS.R")
source ("syntheticData.R")
source ("demons.R")

# Load data from hardrive in .rwl, output it as a TRA ####
load_rwl <- function (stamp=0,fname, ...){
  rwl <- read.rwl(fname, ...)
  
  tra <- rwl.to.tra(rwl)
  tra[tra<=0] <- NA
  
  sdQ <- sample_depth_tra(tra,1)
  sdF <- sample_depth_tra(tra,2)
  sdA <- sample_depth_tra(tra,3)
  
  SD <- list("Q"=sdQ,"F"=sdF,"A"=sdA)
  
  out <- list("tra"=tra,"SD"=SD,"stamp"=stamp)
  
  fname = paste("data/Rdata",stamp,".dat",sep='')
  save(out,file=fname)
  JSONfname = paste("data/JSONdata",stamp,".dat",sep='')
  write_JSON_data(JSONfname,out)
  return (tra) 
}

# Create synthetic data ####
synthetic_tra <- function (stamp=0,method, retentionFraction=1, noiseSD=0, ...){
  if (method=="crude"){
    # Grab the effects from a modern TRA equivalent
    dummy_out <- modern_TRA(...)
    
    effects <- dummy_out$effects
    
    # Rebuild the full TRA
    tra <- effects$I %o% effects$T %o% effects$A
    
    # Keep points at random
    keep <- sample(length(tra), retentionFraction*length(tra))
    tra[!keep] <- NA
    
    # Add noise
    tra <- add_noise (tra, noiseSD, multiplicative)
    
    out <- list(tra=tra, cv=effects)
    
  } else if (method=="modern"){
    out <- modern_TRA(...)
  }
  sdQ <- sample_depth_tra(out$tra,1)
  sdF <- sample_depth_tra(out$tra,2)
  sdA <- sample_depth_tra(out$tra,3)
  
  out$SD <- list("Q"=sdQ,"F"=sdF,"A"=sdA)
  out$stamp <- stamp
  
  fname = paste("data/Rdata",stamp,".dat",sep='')
  save(out,file=fname)
  JSONfname = paste("data/JSONdata",stamp,".dat",sep='')
  write_JSON_data(JSONfname,out)
  
  return (out)
}

# Standardize your tree ring data ####
standardize <- function (dstamp, stamp, ...){
  fname = paste('data/Rdata',dstamp,'.dat',sep='')
  load(fname)
  tra <- out$tra
  out <- standardize_tra(tra, ...)
  nstamp = paste(dstamp,'_',stamp,sep='')
  out$stamp = nstamp
  JSONfname = paste("data/JSONanalysis",dstamp,"_",stamp,".dat",sep='')
  write_JSON_data(JSONfname,out)
  return (out)
  
}

write_JSON_data <- function(fname, out) {
  fileConn<-file(fname)
  writeLines(toJSON(out), fileConn)
  close(fileConn)
  return (0)
}
