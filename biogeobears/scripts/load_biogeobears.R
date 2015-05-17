#!/usr/bin/env Rscript

library(optimx)  # (either 2012 or 2013 version, as of January 2014)
library(FD)        # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; prob. better than library(parallel))
library(BioGeoBEARS)

#######################################################
# SETUP: Load BioGeoBears and Apply Patches
#######################################################

##################################################################################
## This function loads the BioGeoBears package as well as applies various fixes
## and patches required to get it to work correctly.
##
## Based on original made available online by Nick Matzke:
##
##      http://phylo.wikidot.com/biogeobears
##
##
load.biogeobears = function() {
    source("libexec/BioGeoBEARS_add_fossils_randomly_v1.R")
    source("libexec/BioGeoBEARS_basics_v1.R")
    source("libexec/BioGeoBEARS_calc_transition_matrices_v1.R")
    source("libexec/BioGeoBEARS_classes_v1.R")
    source("libexec/BioGeoBEARS_detection_v1.R")
    source("libexec/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
    source("libexec/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
    source("libexec/BioGeoBEARS_generics_v1.R")
    source("libexec/BioGeoBEARS_models_v1.R")
    source("libexec/BioGeoBEARS_on_multiple_trees_v1.R")
    source("libexec/BioGeoBEARS_plots_v1.R")
    source("libexec/BioGeoBEARS_readwrite_v1.R")
    source("libexec/BioGeoBEARS_simulate_v1.R")
    source("libexec/BioGeoBEARS_SSEsim_makePlots_v1.R")
    source("libexec/BioGeoBEARS_SSEsim_v1.R")
    source("libexec/BioGeoBEARS_stochastic_mapping_v1.R")
    source("libexec/BioGeoBEARS_stratified_v1.R")
    source("libexec/BioGeoBEARS_univ_model_v1.R")
    source("libexec/calc_uppass_probs_v1.R")
    source("libexec/calc_loglike_sp_v01.R")
    source("libexec/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
    source("libexec/runBSM_v1.R")
    source("libexec/stochastic_map_given_inputs.R")
    source("libexec/summarize_BSM_tables_v1.R")
    calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
    calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
        # slight speedup hopefully
}
load.biogeobears()

#######################################################
# SETUP: Extension data directory
#######################################################

# When R packages contain extra files, they are stored in the "extdata" directory
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
#
# Each computer operating system might install BioGeoBEARS in a different place,
# depending on your OS and settings.
#
# However, you can find the extdata directory like this:
# extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
# extdata_dir
# list.files(extdata_dir)

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your
# script, since the plot_BioGeoBEARS_results function needs a script from the
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines
# from the package APE which do not pass R CMD check for some reason.
biogeoebears.extdata.dir = normalizePath(system.file("extdata", package="BioGeoBEARS"))

#######################################################
# Convenience: Run Configuration
#######################################################

configure.standard.biogeobears.run = function(bgb.run) {

    # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    bgb.run$speedup=FALSE
    bgb.run$use_optimx = TRUE

    # Multicore processing if desired
    bgb.run$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel"
    # is now default on Macs in R 3.0+, but apparently still
    # has to be typed on some Windows machines. Note: apparently
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)


    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    bgb.run$force_sparse=FALSE

    # Good default settings to get ancestral states
    bgb.run$return_condlikes_table = TRUE
    bgb.run$calc_TTL_loglike_from_condlikes_table = TRUE
    bgb.run$calc_ancprobs = TRUE

    return(bgb.run)
}

get.biogeobears.run = function() {
    bgb.run = define_BioGeoBEARS_run()
    configure.standard.biogeobears.run(bgb.run)
    return(bgb.run)
}

