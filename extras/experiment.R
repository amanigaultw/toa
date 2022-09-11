devtools::install_github("amanigaultw/toa")

rm(list = ls())

library(toa)

#load example data
data("Chang")
data("epith_mesen_ref_raw")

#get diagnosticity scores
toa_ref_epith_mesen <- get_toa_ref(gene_symbols = epith_mesen_ref_raw[,1],
                                   exp_treatment = epith_mesen_ref_raw[,2:11],
                                   exp_control = epith_mesen_ref_raw[,12:21])

#get differentially expressed genes as a function of the stress predictor
DEG_result <- get_DEG(x = Chang$stress,
                      genes = subset(Chang, select = -stress),
                      foldThreshDEG = 1.25)

#get mean diagnosticity scores for these differentially expressed genes
test1 <- toa(x = Chang$stress,
             genes = subset(Chang, select = -stress),
             toa_ref = toa_ref_epith_mesen,
             foldThreshDEG = 1.25)

test1$df.means

x = Chang$stress
genes = subset(Chang, select = -stress)
foldThreshDEG = 1.25
cov = NULL
foldThreshDEG = 1.5
n_boot = 200
show_progress = TRUE

