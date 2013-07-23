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

# Load data from hardrive in .rwl, output it as a TRA ####
load_rwl <- function (stamp=0,fname, ...){
  rwl <- read.rwl(fname, ...)
  
  tra <- rwl.to.tra(rwl)
  
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
synthetic_tra <- function (stamp=0,method, ...){
  if (method=="crude"){
    out <- crude_synth_TRA (...)
    
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
standardize_tra <- function (dstamp, stamp, method, ...){
  fname = paste('data/Rdata',dstamp,'.dat',sep='')
  load(fname)
  tra <- out$tra
  if (method=="RCS"){
    out <- standardize_rcs(tra, ...)
  } else if (method=="SFS"){
    out <- standardize_sfs(tra, ...)#ignore
  } else if (method=="TSFS"){
    out <- standardize_tsfs(tra, ...)
  } else if (method=="FES"){
    out <- standardize_fes(tra, ...)
    out$model <- NULL
  }
  nstamp = paste(dstamp,'_',stamp,sep='')
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

