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
    library(optimx)  # (either 2012 or 2013 version, as of January 2014)
    library(FD)        # for FD::maxent() (make sure this is up-to-date)
    library(snow)     # (if you want to use multicore functionality; prob. better than library(parallel))
    library(BioGeoBEARS)
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
extdata_dir = normalizePath(system.file("extdata", package="BioGeoBEARS"))

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

plot.biogeobears.results.ranges = function(
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

#######################################################
## Convenience: Results Processing (II)
#######################################################

## extracted from wiki

plot.biogeobears.results.areas = function(
        results.object,
        analysis.title="BioGeoBEARS",
        add.corners=FALSE) {

    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=results.object$inputs$geogfn)
    areas = getareas_from_tipranges_object(tipranges)
    max_range_size = results.object$inputs$max_range_size
    trfn = results.object$inputs$trfn
    tr = read.tree(trfn)
    if (is.na(max_range_size)) {
        max_range_size = length(areas)
    }
    states_list_0based = rcpp_areas_list_to_states_list(areas=areas,
            maxareas=max_range_size,
            include_null_range=TRUE)

    barwidth_proportion = 0.02
    barheight_proportion = 0.025
    split_reduce_x = 0.85    # fraction of sizes for the split probabilities
    split_reduce_y = 0.85    # fraction of sizes for the split probabilities

    # Manual offsets, if desired (commented out below)
    offset_nodenums = NULL
    offset_xvals = NULL
    offset_yvals = NULL

    root.edge = TRUE

    # Plot at nodes
    probs_each_area = plot_per_area_probs2(tr,
            res=results.object,
            areas=areas,
            states_list_0based=states_list_0based,
            titletxt=analysis.title,
            cols_each_area=TRUE,
            barwidth_proportion=barwidth_proportion,
            barheight_proportion=barheight_proportion,
            offset_nodenums=offset_nodenums,
            offset_xvals=offset_xvals,
            offset_yvals=offset_yvals,
            root.edge=root.edge)
    # # print("Waiting...")
    # Sys.sleep(2)    # wait for this to finish

    # Calculate per-area probabilities for corners
    relprobs_matrix = results.object$ML_marginal_prob_each_state_at_branch_bottom_below_node
    probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)

    if (add.corners) {
        # Add to left corners
        #offset_nodenums = c(20, 21)
        #offset_xvals = c(0, 0)
        #offset_yvals = c(-1, -0.25)
        add_per_area_probs_to_corners2(tr,
                areas,
                probs_each_area,
                left_or_right="left",
                cols_each_area=TRUE,
                barwidth_proportion=split_reduce_x*barwidth_proportion,
                barheight_proportion=split_reduce_y*barheight_proportion,
                offset_nodenums=offset_nodenums,
                offset_xvals=offset_xvals,
                offset_yvals=offset_yvals,
                root.edge=root.edge)
        # print("Waiting...")
    }
}

#### BELOW OVERRIDES BIOGEOBEARS LIBRARY CODE ####

plot_per_area_probs2 = function(tr, res, areas, states_list_0based, titletxt="", cols_each_area=NULL, barwidth_proportion=0.02, barheight_proportion=0.025, offset_nodenums=NULL, offset_xvals=NULL, offset_yvals=NULL, root.edge=TRUE, border="default", trcol="black", plot_rangesizes=FALSE) {
	defaults='
	areas = getareas_from_tipranges_object(tipranges)
	states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)
	states_list_0based

	cols_each_area=NULL
	offset_nodenums=NULL
	offset_xvals=NULL
	offset_yvals=NULL
	barwidth_proportion=0.02
	barheight_proportion=0.025

	border="default"
	trcol="black"
	plot_rangesizes=FALSE
	' # END defaults


	# Error check
	if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
		{
		errortxt = "\n\nERROR IN plot_per_area_probs2():\ncols_each_area=FALSE, but it should be either:\nNULL (which gives grey boxes)\nTRUE (which makes colors auto-generated), or \na list of colors the same length as 'areas'.\n\n"
		cat(errortxt)

		stop("\n\nStopping on error.\n\n")

		} # END if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))


	# Get the relative probabilities of each state/range
	relprobs_matrix = res$ML_marginal_prob_each_state_at_branch_top_AT_node

	# Collapse to the probabilities of each area
	if (plot_rangesizes == FALSE)
		{
		probs_each_area = infprobs_to_probs_of_each_area(relprobs_matrix, states_list=states_list_0based)
		}

	# Collapse to the probabilities of each range size
	if (plot_rangesizes == TRUE)
		{
		rangesizes_by_state = sapply(FUN=length, X=states_list_0based)
		rangesizes = sort(unique(rangesizes_by_state))
		rangesizes

		probs_each_area = infprobs_to_probs_of_each_rangesize(relprobs_matrix, states_list=states_list_0based)
		}
	dim(probs_each_area)



	# To get offset_tiplabels:
	ntips = length(tr$tip.label)
	numnodes = tr$Nnode
	tipnums = 1:ntips
	nodenums = (ntips+1):(ntips+numnodes)


	nodes_xy = node_coords(tr,
	                       root.edge=root.edge,
                           tmplocation=paste(extdata_dir, "a_scripts" , sep="/"))
	nodes_xy

	# Make bar width a proportion of the width of the plot in x
	#barwidth_proportion = 0.02
	barwidth = (max(nodes_xy$x) - min(nodes_xy$x)) * barwidth_proportion
	barwidth

	#barheight_proportion = 0.025
	barheight = (max(nodes_xy$y) - min(nodes_xy$y)) * barheight_proportion
	barheight

	numareas = ncol(probs_each_area)
	areanums = 1:numareas
	middle = median(areanums)
	middle
	offsets_nodes = (areanums - middle)*barwidth
	offsets_tips = (areanums - 0)*barwidth
	offset_tiplabels = max(offsets_tips) + barwidth/1

	# Plot the tree
	plot(tr, label.offset=offset_tiplabels, root.edge=root.edge, edge.color=trcol)
	# plot(tr, label.offset=offset_tiplabels)
	axisPhylo()
	title(titletxt)


	# Add the areas boxes
	add_per_area_probs_to_nodes2(tr,
	        areas=areas,
            probs_each_area,
            cols_each_area=cols_each_area,
            barwidth_proportion=barwidth_proportion,
            barheight_proportion=barheight_proportion,
            offset_nodenums=offset_nodenums,
            offset_xvals=offset_xvals,
            offset_yvals=offset_yvals,
            border=border)

	return(probs_each_area)
	}


add_per_area_probs_to_corners2 = function(tr, areas, probs_each_area, left_or_right, cols_each_area=NULL, barwidth_proportion=0.02, barheight_proportion=0.025, offset_nodenums=NULL, offset_xvals=NULL, offset_yvals=NULL, root.edge=TRUE, border="default", trcol)
	{
	defaults='
	cols_each_area=NULL
	offset_nodenums=NULL
	offset_xvals=NULL
	offset_yvals=NULL
	barwidth_proportion=0.02
	barheight_proportion=0.025
	' # END defaults

	# Error check
	if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
		{
		errortxt = "\n\nERROR IN add_per_area_probs_to_corners2():\n\ncols_each_area=FALSE, but it should be either:\nNULL (which gives grey boxes)\nTRUE (which makes colors auto-generated), or \na list of colors the same length as 'areas'.\n\n"
		cat(errortxt)

		stop("\n\nStopping on error.\n\n")

		} # END if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))


	ntips = length(tr$tip.label)
	numnodes = tr$Nnode
	tipnums = 1:ntips
	#nodenums = (ntips+1):(ntips+numnodes)
	nodenums = (ntips+1):(ntips+numnodes)

	# Get the plot coordinates of the corners ABOVE each node
	# corners_list = corner_coords(tr, root.edge=root.edge)
	corners_list = corner_coords(tr,
	                       root.edge=root.edge,
                           tmplocation=paste(extdata_dir, "a_scripts" , sep="/"))
	corners_list

	# Get the node numbers of the nodes ABOVE each corner above each node
	leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr)
	leftright_nodes_matrix

	if (left_or_right == "left")
		{
		corners_xy = corners_list$leftcorns
		nodenums_above_LorR_corner = leftright_nodes_matrix$left
		}

	if (left_or_right == "right")
		{
		corners_xy = corners_list$rightcorns
		nodenums_above_LorR_corner = leftright_nodes_matrix$right
		}




	# LEFT OR RIGHT SPLITS
	# Probs of each area BELOW the node (nodenum = rownum)

	# Make bar width a proportion of the width of the plot in x
	#barwidth_proportion = 0.02
	barwidth = (max(corners_xy$x) - min(corners_xy$x)) * barwidth_proportion
	barwidth

	#barheight_proportion = 0.025
	barheight = (max(corners_xy$y) - min(corners_xy$y)) * barheight_proportion
	barheight

	numareas = ncol(probs_each_area)
	areanums = 1:numareas
	middle = median(areanums)
	middle
	offsets_nodes = (areanums - middle)*barwidth
	offsets_tips = (areanums - 0)*barwidth
	offset_tiplabels = max(offsets_tips) + barwidth/1

	# Draw the empty boxes

	# xcoords for tips
	nodenums_above_LorR_corner

	# We just need to plot at the corners above internal nodes, not tips
	#xlefts_tips = sapply(X=corners_xy$x[nodenums], FUN="+", (offsets_tips - barwidth/2))
	rownums = nodenums - ntips
	xlefts_nodes = sapply(X=corners_xy$x[rownums], FUN="+", (offsets_nodes - barwidth/2))
	#xrights_tips = sapply(X=corners_xy$x[tipnums], FUN="+", (offsets_tips + barwidth/2))
	xrights_nodes = sapply(X=corners_xy$x[rownums], FUN="+", (offsets_nodes + barwidth/2))

	xlefts_matrix = t(xlefts_nodes)
	xrights_matrix = t(xrights_nodes)

	ybottoms_per_node = sapply(X=corners_xy$y, FUN="-", (barheight/2))
	ytops_per_node = sapply(X=corners_xy$y, FUN="+", (barheight/2))

	# Manually modify some positions
	if (is.null(offset_nodenums) == FALSE)
		{
		rownums = match(x=offset_nodenums, table=nodenums)
		xlefts_matrix[rownums,] = xlefts_matrix[rownums,] + (offset_xvals*barwidth)
		xrights_matrix[rownums,] = xrights_matrix[rownums,] + (offset_xvals*barwidth)
		ybottoms_per_node[rownums] = ybottoms_per_node[rownums] + (offset_yvals*barheight)
		ytops_per_node[rownums] = ytops_per_node[rownums] + (offset_yvals*barheight)
		}

	# Convert these into plain lists
	xlefts = c(xlefts_matrix)
	xrights = c(xrights_matrix)

	ybottoms = rep(ybottoms_per_node, times=numareas)
	ytops = rep(ytops_per_node, times=numareas)


	# Plot the box outlines
	#nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

	# Then draw black boxes inside these
	# You just have to adjust the top of the black bar, based on prob
	yranges_per_node = ytops_per_node - ybottoms_per_node
	yadd_above_ybottom = yranges_per_node * probs_each_area[nodenums_above_LorR_corner,]
	ytops_probs_node = yadd_above_ybottom + ybottoms
	ytops_probs = c(ytops_probs_node)


	# IF cols_each_area == TRUE, auto-generate colors
	if (length(cols_each_area) == 1 && (cols_each_area == TRUE))
		{
		tmp_colors_matrix = get_colors_for_numareas(numareas, use_rainbow=FALSE)
		#cols_each_area = c("blue", "green", "yellow", "red")
		cols_each_area = mapply(FUN=rgb, red=tmp_colors_matrix[1,], green=tmp_colors_matrix[2,], blue=tmp_colors_matrix[3,], MoreArgs=list(maxColorValue=255))
		}


	# Default color is darkgray
	if (is.null(cols_each_area))
		{
		# Auto-generate border, also -- gray50 looks black against lighter gray70
		if (border == "default")
			{
			border = "gray50"
			}

		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

		# Draw bars
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, MoreArgs=list(col="gray70", border=border))
		} else {


		if ( length(cols_each_area) != length(areas))
			{
			errortxt = paste("\n\nERROR in add_per_area_probs_to_corners2():\n\nif cols_each_area is specified, length(cols_each_area) must equal length(areas).\n\n", sep="")
			cat(errortxt)

			cat("Your 'areas':\n\n", sep="")
			print(areas)

			cat("Your 'cols_each_area':\n\n", sep="")
			print(cols_each_area)

			stop("Stopping on error.")
			} # END if ( length(cols_each_area) != length(area))

		# Otherwise, expand colors to each box
		cols_each_area_matrix = matrix(data=cols_each_area, nrow=numnodes, ncol=length(cols_each_area), byrow=TRUE)


		# Plot with different internal box colors
		cols_each_area = c(cols_each_area_matrix)

		# Auto-generate border with colored boxes (black looks best), also
		if (border == "default")
			{
			border = "black"
			}

		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, col=cols_each_area, MoreArgs=list(border=border))
		} # END if (is.null(cols_each_area))

	return(NULL)
	}


add_per_area_probs_to_nodes2 = function(tr, areas, probs_each_area, cols_each_area=NULL, barwidth_proportion=0.02, barheight_proportion=0.025, offset_nodenums=NULL, offset_xvals=NULL, offset_yvals=NULL, root.edge=TRUE, border="default")
	{
	defaults='
	cols_each_area=NULL
	offset_nodenums=NULL
	offset_xvals=NULL
	offset_yvals=NULL
	barwidth_proportion=0.02
	barheight_proportion=0.025
	' # END defaults


	# Error check
	if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))
		{
		errortxt = "\n\nERROR IN add_per_area_probs_to_nodes2():\n\ncols_each_area=FALSE, but it should be either:\nNULL (which gives grey boxes)\nTRUE (which makes colors auto-generated), or \na list of colors the same length as 'areas'.\n\n"
		cat(errortxt)

		stop("\n\nStopping on error.\n\n")

		} # END if ((length(cols_each_area) == 1) && (cols_each_area == FALSE))


	ntips = length(tr$tip.label)
	numnodes = tr$Nnode
	tipnums = 1:ntips
	nodenums = (ntips+1):(ntips+numnodes)

	# nodes_xy = node_coords(tr, root.edge=root.edge)
	nodes_xy = node_coords(tr,
	                       root.edge=root.edge,
                           tmplocation=paste(extdata_dir, "a_scripts" , sep="/"))
	nodes_xy

	# Make bar width a proportion of the width of the plot in x
	#barwidth_proportion = 0.02
	barwidth = (max(nodes_xy$x) - min(nodes_xy$x)) * barwidth_proportion
	barwidth

	#barheight_proportion = 0.025
	barheight = (max(nodes_xy$y) - min(nodes_xy$y)) * barheight_proportion
	barheight

	numareas = ncol(probs_each_area)
	areanums = 1:numareas
	middle = median(areanums)
	middle
	offsets_nodes = (areanums - middle)*barwidth
	offsets_tips = (areanums - 0)*barwidth
	offset_tiplabels = max(offsets_tips) + barwidth/1

	# Draw the empty boxes

	# xcoords for tips
	xlefts_tips = sapply(X=nodes_xy$x[tipnums], FUN="+", (offsets_tips - barwidth/2))
	xlefts_nodes = sapply(X=nodes_xy$x[nodenums], FUN="+", (offsets_nodes - barwidth/2))
	xrights_tips = sapply(X=nodes_xy$x[tipnums], FUN="+", (offsets_tips + barwidth/2))
	xrights_nodes = sapply(X=nodes_xy$x[nodenums], FUN="+", (offsets_nodes + barwidth/2))

	xlefts_matrix = t(cbind(xlefts_tips, xlefts_nodes))
	xrights_matrix = t(cbind(xrights_tips, xrights_nodes))

	ybottoms_per_node = sapply(X=nodes_xy$y, FUN="-", (barheight/2))
	ytops_per_node = sapply(X=nodes_xy$y, FUN="+", (barheight/2))

	# Manually modify some positions
	if (is.null(offset_nodenums) == FALSE)
		{
		xlefts_matrix[offset_nodenums,] = xlefts_matrix[offset_nodenums,] + (offset_xvals*barwidth)
		xrights_matrix[offset_nodenums,] = xrights_matrix[offset_nodenums,] + (offset_xvals*barwidth)
		ybottoms_per_node[offset_nodenums] = ybottoms_per_node[offset_nodenums] + (offset_yvals*barheight)
		ytops_per_node[offset_nodenums] = ytops_per_node[offset_nodenums] + (offset_yvals*barheight)
		}

	# Convert these into plain lists
	xlefts = c(xlefts_matrix)
	xrights = c(xrights_matrix)

	ybottoms = rep(ybottoms_per_node, times=numareas)
	ytops = rep(ytops_per_node, times=numareas)


	# Plot the box outlines
	#nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

	# Then draw black boxes inside these
	# You just have to adjust the top of the black bar, based on prob
	yranges_per_node = ytops_per_node - ybottoms_per_node
	yadd_above_ybottom = yranges_per_node * probs_each_area
	ytops_probs_node = yadd_above_ybottom + ybottoms
	ytops_probs = c(ytops_probs_node)


	# IF cols_each_area == TRUE, auto-generate colors
	if (length(cols_each_area) == 1 && (cols_each_area == TRUE))
		{
		tmp_colors_matrix = get_colors_for_numareas(numareas, use_rainbow=FALSE)
		#cols_each_area = c("blue", "green", "yellow", "red")
		cols_each_area = mapply(FUN=rgb, red=tmp_colors_matrix[1,], green=tmp_colors_matrix[2,], blue=tmp_colors_matrix[3,], MoreArgs=list(maxColorValue=255))

		}



	# Default color is darkgray
	if (is.null(cols_each_area))
		{
		# Auto-generate border, also -- gray50 looks black against lighter gray70
		if (border == "default")
			{
			border = "gray50"
			}

		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

		# Draw bars
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, MoreArgs=list(col="gray70", border=border))
		} else {


		if ( length(cols_each_area) != length(areas))
			{
			errortxt = paste("\n\nERROR in add_per_area_probs_to_nodes2():\n\nif cols_each_area is specified, length(cols_each_area) must equal length(areas).\n\n", sep="")

			cat(errortxt)

			cat("Your 'areas':\n\n", sep="")
			print(areas)

			cat("Your 'cols_each_area':\n\n", sep="")
			print(cols_each_area)

			stop("Stopping on error.")
			} # END if ( length(cols_each_area) != length(areas))

		# Otherwise, expand colors to each box
		cols_each_area_matrix = matrix(data=cols_each_area, nrow=(ntips+numnodes), ncol=length(cols_each_area), byrow=TRUE)


		# Plot with different internal box colors
		cols_each_area = c(cols_each_area_matrix)

		# Auto-generate border with colored boxes (black looks best), also
		if (border == "default")
			{
			border = "black"
			}

		# Plot the box outlines
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops, MoreArgs=list(col="white", border=border))

		# Plot the colored bars
		nulls = mapply(FUN=rect, xleft=xlefts, ybottom=ybottoms, xright=xrights, ytop=ytops_probs, col=cols_each_area, MoreArgs=list(border=border))
		} # END if (is.null(cols_each_area))

	return(NULL)
	}
