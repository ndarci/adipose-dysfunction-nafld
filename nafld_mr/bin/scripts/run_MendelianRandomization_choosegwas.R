library(MendelianRandomization)

# import input data in TwoSampleMR format
mr = read.table('../../data/mr_input_table_TwoSampleMR.txt')

# convert to MendelianRandomization format
MRInputObject = mr_input(bx = mr$beta.exposure, 
                          bxse = mr$se.exposure,
                          by = mr$beta.outcome, 
                          byse = mr$se.outcome)

# set up dataframe to store important stats from each method
res.cols = c('method', 'beta', 'p.mr', 'stat.hetero', 'p.hetero')
res = data.frame(matrix(NA, nrow = 4, ncol = length(res.cols)))
names(res) = res.cols

# run IVW MR method
IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)
print(IVWObject)
res[1,] = c('IVW', IVWObject$Estimate, IVWObject$Pvalue, 
            IVWObject$Heter.Stat[1], IVWObject$Heter.Stat[2])

# also run IVW for just the VEGFB cis regulatory IV
mr.vegfb = mr[mr$SNP == 'rs56271783',]
MRInputObject.vegfb = mr_input(bx = mr.vegfb$beta.exposure, 
                          bxse = mr.vegfb$se.exposure,
                          by = mr.vegfb$beta.outcome, 
                          byse = mr.vegfb$se.outcome)
IVWObject.vegfb <- mr_ivw(MRInputObject.vegfb,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = FALSE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05)
res[2,] = c('IVW_VEGFB_only', IVWObject.vegfb$Estimate, IVWObject.vegfb$Pvalue, 
            NA, NA)

# run MR-Egger method
EggerObject <- mr_egger(MRInputObject, 
                        robust = FALSE,
                        penalized = FALSE,
                        correl = FALSE,
                        distribution = "normal",
                        alpha = 0.05)
print(EggerObject)
res[3,] = c('MR-Egger', EggerObject$Estimate, EggerObject$Pvalue.Est,
            EggerObject$Heter.Stat[1], EggerObject$Heter.Stat[2])

# run weighted median method
WeightedMedianObject <- mr_median(MRInputObject, 
                                  weighting = "weighted", 
                                  distribution = "normal", 
                                  alpha = 0.05, 
                                  iterations = 10000, 
                                  seed = 314159265)
print(WeightedMedianObject)
res[4,] = c('Weighted_Median', WeightedMedianObject$Estimate, WeightedMedianObject$Pvalue,
            NA, NA)

write.table(res, '../../data/other_mr_methods_results_choosegwas.txt')











