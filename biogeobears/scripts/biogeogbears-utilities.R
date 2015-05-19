#!/usr/bin/env Rscript

###############################################################################
##
##  Copyright 2015 Jeet Sukumaran.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
##  NOTE: Some of the following code has been based, derived, modified
##  from or otherwise taken entirely from code released by Nick Matzke
##  and made available (under the GPL3 license) at:
##
##      http://phylo.wikidot.com/biogeobears
##
##  and are:
##
##      Copyright (C) Nick Matzke.
##
##
###############################################################################

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
## Convenience: Run Configuration
##
## Based on original made available online by Nick Matzke:
##
##      http://phylo.wikidot.com/biogeobears
##
#######################################################

configure.standard.biogeobears.run = function(results.object) {

    # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    results.object$speedup=FALSE
    results.object$use_optimx = TRUE

    # Multicore processing if desired
    results.object$num_cores_to_use = 1
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
    results.object$force_sparse=FALSE

    # Good default settings to get ancestral states
    results.object$return_condlikes_table = TRUE
    results.object$calc_TTL_loglike_from_condlikes_table = TRUE
    results.object$calc_ancprobs = TRUE

    return(results.object)
}

get.biogeobears.run = function() {
    results.object = define_BioGeoBEARS_run()
    configure.standard.biogeobears.run(results.object)
    return(results.object)
}

#######################################################
## Convenience: Results Processing
#######################################################

plot.biogeobears.results = function(
        results.object,
        plot.type="text",
        analysis.title="BioGeoBEARS") {
    res2 = plot_BioGeoBEARS_results(
            results.object,
            analysis.title,
            addl_params=list("j"),
            plotwhat=plot.type,
            label.offset=0.45,
            tipcex=0.7,
            statecex=0.7,
            splitcex=0.6,
            titlecex=0.8,
            plotsplits=TRUE,
            cornercoords_loc=normalizePath(system.file("extdata/a_scripts", package="BioGeoBEARS")),
            include_null_range=TRUE,
            tr=read.tree(results.object$inputs$trfn),
            tipranges=getranges_from_LagrangePHYLIP(lgdata_fn=results.object$inputs$geogfn)
            )
    res2
}

get.biogeobears.ranges = function(geogfn, max.range.size=NULL) {
    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=normalizePath(geogfn))
    areas = getareas_from_tipranges_object(tipranges)
    if (is.null(max.range.size) || is.na(max.range.size)) {
        max.range.size = length(areas)
    }
    state_indices_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max.range.size, include_null_range=TRUE)
    ranges_list = areas_list_to_states_list_new(areas=areas, maxareas=length(areas), include_null_range=TRUE, split_ABC=FALSE)
    ranges = unlist(ranges_list)
    rval = list(
                tipranges=tipranges,
                areas=areas,
                state_indices=state_indices_0based,
                ranges_list=ranges_list,
                ranges=ranges)
    return(rval)
}

get.biogeobears.results.table = function(results.object) {
    # Get the likelihood of the ranges as a data.frame
    # This has the likelihood of the states (ranges) of each node as columns,
    # with each row listing the node whose index corresponds to the row index
    results1 = data.frame(dec.results$ML_marginal_prob_each_state_at_branch_top_AT_node)

    # Rename the columns with range names
    geogfn = results.object$inputs$geogfn
    max.range.size = results.object$inputs$max_range_size
    geog.info = get.biogeobears.ranges(geogfn, max.range.size=max.range.size)
    names(results1) <-  geog.info$ranges

    # Add the node info
    trfn = results.object$inputs$trfn
    tree = read.tree(trfn)
    results2 = cbind(prt(tree, get_tipnames=TRUE), results1)
    results2$daughter_nds = NULL # elements are a list (might be a better way to handle this?)
    return(results2)
}
