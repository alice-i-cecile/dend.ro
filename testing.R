# Trying out existing data sets ####
library(dplR)
data(ca533); rwl <- ca533
data(co021); rwl <- co021
data(gp.rwl); rwl <- gp.rwl
data(gp.po); rwl <- bai.in(gp.rwl,gp.po)

# Remove 0s
rwl[rwl==0] <- NA

tra <- rwl.to.tra(rwl)

std.RCS.AT <- standardize_tra(tra, model=list(I=F, A=T, T=T), method="rcs")
fit.RCS.AT <- unlist(std.RCS.AT$fit[3:13])
rm(std.RCS.AT)

std.RCS.TA <- standardize_tra(tra, model=list(I=F, T=T, A=T), method="rcs")
fit.RCS.TA <- unlist(std.RCS.TA$fit[3:13])
rm(std.RCS.TA)

std.RCS.ITA <- standardize_tra(tra, model=list(I=T, T=T, A=T), method="rcs")
fit.RCS.ITA <- unlist(std.RCS.ITA$fit[3:13])
rm(std.RCS.ITA)

std.SFS.AT <- standardize_tra(tra, model=list(I=F, A=T, T=T), method="sfs")
fit.SFS.AT <- unlist(std.SFS.AT$fit[3:13])
rm(std.SFS.AT)

std.SFS.TA <- standardize_tra(tra, model=list(I=F, T=T, A=T), method="sfs")
fit.SFS.TA <- unlist(std.SFS.TA$fit[3:13])
rm(std.SFS.TA)

std.SFS.ITA <- standardize_tra(tra, model=list(I=T, T=T, A=T), method="sfs")
fit.SFS.ITA <- unlist(std.SFS.ITA$fit[3:13])
rm(std.SFS.ITA)

fit_df <- data.frame(
  fit.RCS.AT,
  fit.RCS.TA,
  fit.RCS.ITA,
  fit.SFS.AT,
  fit.SFS.TA,
  fit.SFS.ITA
)

print(fit_df)

# AIC Weights
aic_weights <- function(AIC)
{
  dAIC <- AIC - min(AIC)
  weights <- exp(-0.5*dAIC)/sum(exp(-0.5*dAIC))
}

print(aic_weights(fit_df["AIC",]))
print(aic_weights(fit_df["AICc",]))
print(aic_weights(fit_df["BIC",]))

# Format conversion ####
tra <- rwl.to.tra(rwl)
stra <- sparse_tra(tra)
tra2 <- unsparse_tra(stra)
rwl2 <- tra.to.rwl(tra2)
stra2 <- rwl.to.stra(rwl)
tra3 <- unsparse_tra(stra2)

# Performance
system.time(rwl.to.tra(rwl))
system.time(sparse_tra(rwl.to.tra(rwl)))
system.time(rwl.to.stra(rwl))

# Sparse testing ####
full_rcs <- standardize_tra(tra, sparse=F, method="rcs")
sp_rcs <- standardize_tra(tra, sparse=T, method="rcs")
full_sfs <- standardize_tra(tra, sparse=F, method="sfs")
sp_sfs <- standardize_tra(tra, sparse=T, method="sfs")

system.time(standardize_tra(tra, sparse=F, method="rcs"))
system.time(standardize_tra(tra, sparse=T, method="rcs"))
system.time(standardize_tra(tra, sparse=F, method="sfs"))
system.time(standardize_tra(tra, sparse=T, method="sfs"))


fit_full_rcs <- unlist(full_rcs$fit[3:13])
fit_sp_rcs <- unlist(sp_rcs$fit[3:13])
fit_full_sfs <- unlist(full_sfs$fit[3:13])
fit_sp_sfs <- unlist(sp_sfs$fit[3:13])

fit_df <- data.frame(
  fit_full_rcs,
  fit_sp_rcs,
  fit_full_sfs,
  fit_sp_sfs
)

print(fit_df)