#!/usr/bin/env Rscript

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
# We have to extract the log-likelihood differently, depending on the
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(dec.results)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(decj.results)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(
        results_object=dec.results,
        returnwhat="table",
        addl_params=c("j"),
        paramsstr_digits=4)
# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(
        results_object=decj.results,
        returnwhat="table",
        addl_params=c("j"),
        paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
# confer the same likelihood on the data. See: Brian O'Meara's webpage:
# http://www.brianomeara.info/tutorials/aic
# ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

teststable$alt = c("DEC+J")
teststable$null = c("DEC")
row.names(restable) = c("DEC", "DEC+J")
