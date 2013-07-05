linear_trend <- function(x, k=1, s=0.01){return (k+s*x)}

test <- modern_TRA(nQ=10, nF=10, sdF=0.5, funcA=linear_trend, noiseSD=0.5, GSS=make_GSS(s=0.5, lambda=50, k=50, Beta=0))

lm_test <- standardize_fes (test$tra, incF=F, incA=F, model_type="lm", multiplicative=T)

gnm_test <- standardize_fes (test$tra, incA=F, model_type="gnm", multiplicative=T)

gam_test <- standardize_fes (test$tra, incA=F, model_type="gam", multiplicative=T)

rcs_test <- standardize_rcs(test$tra, factor_order=c(3), meanType="arithmetic")

rsfs_test <- standardize_sfs(test$tra, factor_order=c(3,1, 2))

tsfs_test <- standardize_tsfs(test$tra, factor_order=c(3,1, 2))


standardize_rcs(test$tra, factor_order=c(3), meanType="arithmetic")
standardize_rcs(test$tra, factor_order=c(3), meanType="geometric")