
#######################################################
# Run BioGeoBEARS on multiple trees
#######################################################

#' @BioGeoBEARS_run_object Set up the inference model you want to run on 
#' each tree, like you would for a normal single-tree run.
#' @newick_fns A list of the Newick files (e.g., you should extract some trees
#' from a BEAST NEXUS MCMC output)
#' @model_name The name you would like added to output filenames
#' @geog_fns A list of corresponding geography files (by default, these are just .geog instead of .newick)
#' @resfns A list of results filenames for .Rdata files, either for saving BioGeoBEARS analyses on each tree, or 
#' for loading previously-saved runs (by default, these are just _model_name.Rdata instead of .newick)
#' @run_or_collect If you want to run BioGeoBEARS on each tree (slow), use "run". If you just want to 
#' collect the results over all the trees, use "collect".  For both, pick "both".
#' @start_treenum Default 1. Change if you want to skip some trees
#' @end_treenum Default is length(newick_fns). Change if you want to run a subset of trees.
#' @runslow If FALSE, old, saved .Rdata results are loaded via load(). Default TRUE generates new .Rdata files and 
#' saves them via save()
#' 

run_bears_optim_on_multiple_trees <- function(BioGeoBEARS_run_object, newick_fns, model_name="", geog_fns=NULL, resfns=NULL, run_or_collect="collect", start_treenum=1, end_treenum=length(newick_fns), runslow=TRUE, plot_params=FALSE)
	{
	# Check the input newick_fns to make sure they end in .newick
	num_newick_strings = str_count(string=newick_fns, pattern="\\.newick")
	num_newick_strings_equals_1_TF = num_newick_strings == 1
	if (sum(num_newick_strings_equals_1_TF) != length(num_newick_strings_equals_1_TF))
		{
		error_txt = paste("\n\nERROR in run_bears_optim_on_multiple_trees(): All filenames in 'newick_fns' must have one and only one '.newick'.\nViolators in your 'newick_fns':\n\n", sep="")
		cat(error_txt)
		cat(newick_fns[num_newick_strings_equals_1_TF=FALSE], sep="\n")
		
		stop("Stopping on error in run_bears_optim_on_multiple_trees()")
		}
	
	# File names to store the state probabilities at nodes and corners
	if (model_name == "")
		{
		suffix_txt = ""
		} else {
		suffix_txt = paste("_", model_name, sep="")
		}
	stateprobs_at_nodes_fn = paste("state_probs_at_nodes_across_all_trees", suffix_txt, ".Rdata", sep="")
	stateprobs_at_corners_fn = paste("state_probs_at_corners_across_all_trees", suffix_txt, ".Rdata", sep="")
	
	# Get names of:
	# - geography filenames (geog_fns)
	# - results filenames (resfns)
	#newick_fns = output_trfns
	if (is.null(geog_fns))
		{
		geog_fns = gsub(pattern="newick", replacement="geog", x=newick_fns)
		cat(geog_fns, ",")
		}
	if (is.null(resfns))
		{
		replacement_txt = paste(suffix_txt, ".Rdata", sep="")
		resfns = gsub(pattern=".newick", replacement=replacement_txt, x=newick_fns)
		}

	if (runslow == TRUE)
		{
		row_start = 1
		row_end = 0
		for (i in start_treenum:end_treenum)
		#for (i in 1:length(newick_fns))
		#for (i in 1:2)
			{
			# If you just want to run the inferences and save the results
			if ((run_or_collect == "run") || (run_or_collect == "both"))
				{
				txt = paste("\n\nRunning inference under '", model_name, "' for tree #", i, " of ", start_treenum, "-", end_treenum, "...\n\n", sep="")
				cat(txt)
		
				BioGeoBEARS_run_object$geogfn = geog_fns[i]
				BioGeoBEARS_run_object$trfn = newick_fns[i]
				resfn = resfns[i]
		
				# Check the max number of areas:
				#tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn)
				# max(rowSums(dfnums_to_numeric(tipranges@df), na.rm=TRUE))

				res = bears_optim_run(BioGeoBEARS_run_object)
				res    
		
				save(res, file=resfn)
				} # end if (run_or_collect == TRUE)
			
			# Initialize the matrices if it's the first iteration
			if (i == start_treenum)
				{
				# If collecting results, create big empty tables to store the state probabilities at nodes and corners
				if ((run_or_collect == "collect") || (run_or_collect == "both"))
					{
					# For summarizing over the tree collection:
					# Set up an empty table to hold the state
					# probabilities of each run
					# numrows = numtrees * numnodes per tree
					# numcols = numstates
					example_trfn = newick_fns[[start_treenum]]
					example_tr = read.tree(example_trfn)
					numrows = length(example_tr$tip.label) + example_tr$Nnode
					numrows
					# Load results to res
					load(resfns[[1]])
					numstates = ncol(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
	
					# We will make these matrices double-size, just in case trees vary in size
					# (e.g., fossil inclusion or not), then cut the matrices down by excluding 
					# all NAs in the last column.
	
					# Create a big empty matrix (+1 to numcols for the OTUnames)
					# for the node states
					state_probs_at_nodes_across_all_trees = matrix(data=NA, nrow=length(resfns)*numrows*2, ncol=numstates+1)
					dim(state_probs_at_nodes_across_all_trees)
		
					# Create another matrix for the states at the corners (each corner is below
					# a node; thus probably nothing below the root)
					numrows = nrow(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
					numstates = ncol(res$ML_marginal_prob_each_state_at_branch_bottom_below_node)
					state_probs_at_corners_across_all_trees = matrix(data=NA, nrow=length(resfns)*numrows*2, ncol=numstates+1)
					dim(state_probs_at_corners_across_all_trees)
					} # end if ((run_or_collect == "collect") || (run_or_collect == "both"))
				} # end if (i == start_treenum)
		

	
	
			# If you want to run the summary after each tree, or just the summaries
			if ((run_or_collect == "collect") || (run_or_collect == "both"))
				{
				# Processing previously done inferences
				txt = paste("\nProcessing previous inferences under '", model_name, "' for tree #", i, " of ", start_treenum, "-", end_treenum, "...", sep="")
				cat(txt)


				# Do processing if desired
				# The goal is:
				# For each node and corner, get:
				# 1. The number of times that node appears in the sample
				#    (this should be approximately the PP of the bipartition)
				# 2. Whenever the node does appear, get the state probabilities
				#    at that node.
				# 3. Sum over all state probabilities and divide by #1
				# 4. Result is ancestral state probabilities averaged over the tree
		
				# For each tree, make an sorted text list of the OTUs descending from each node
				# Add to the state probabilities for each node
				trfn = newick_fns[i]
				tr = read.tree(trfn)
				trtable = prt(tr, printflag=FALSE, get_tipnames=TRUE)
				head(trtable)

				# Get the BioGeoBEARS results for this tree
				# Loads to "res"
				resfn = resfns[i]
				load(resfn)
		
				names(res)
				state_probs_at_nodes = res$ML_marginal_prob_each_state_at_branch_top_AT_node
				state_probs_at_corners = res$ML_marginal_prob_each_state_at_branch_bottom_below_node
				dim(state_probs_at_corners)
				dim(trtable)
		
				# Store the results
				row_end = row_start-1 + nrow(state_probs_at_nodes)
				rownums = row_start:row_end
				# first columns are the state probabilities
				colnums = 1:numstates
				names_col = numstates + 1
		
				# Store the state probabilities at nodes
				state_probs_at_nodes_across_all_trees[rownums, colnums] = state_probs_at_nodes
				state_probs_at_corners_across_all_trees[rownums, colnums] = state_probs_at_corners
		
				# Store the OTUs descending from each node
				OTUnames = trtable$tipnames
				state_probs_at_nodes_across_all_trees[rownums, names_col] = OTUnames
				state_probs_at_corners_across_all_trees[rownums, names_col] = OTUnames
				
				# update row_start
				row_start = row_end + 1
				} # end if ((run_or_collect == FALSE) or (run_or_collect == "both"))
			} # end for (i in 1:length(newick_fns))
		
		# Delete rows that are all NA
		allNA <- function(tmprow)
			{
			row_is_allNA_TF = all(is.na(tmprow))
			return(row_is_allNA_TF)
			}
		rows_allNA_TF = apply(X=state_probs_at_nodes_across_all_trees, MARGIN=1, FUN=allNA)
		state_probs_at_nodes_across_all_trees = state_probs_at_nodes_across_all_trees[rows_allNA_TF==FALSE,]
		state_probs_at_corners_across_all_trees = state_probs_at_corners_across_all_trees[rows_allNA_TF==FALSE,]
		
		# Save the states		
		txt = paste("\n\nSaving state probabilities at nodes and corners across your tree sample to:\nworking directory: ", getwd(), "\n'", stateprobs_at_nodes_fn, "'\n'", stateprobs_at_corners_fn, "'\n(may be slow, ~1 minute)\n", sep="")
		cat(txt)
		
		# After the for-loop, save the ancestral states matrices if you like
		save(state_probs_at_nodes_across_all_trees, file=stateprobs_at_nodes_fn)
		save(state_probs_at_corners_across_all_trees, file=stateprobs_at_corners_fn)
		} else {
		
		# Load the states
		txt = paste("\n\nLoading state probabilities at nodes and corners across your tree sample from:\nworking directory: ", getwd(), "\n'", stateprobs_at_nodes_fn, "'\n'", stateprobs_at_corners_fn, "'\n(may be slow, ~1 minute)\n", sep="")
		cat(txt)
		# Load to: state_probs_at_nodes_across_all_trees
		load(file=stateprobs_at_nodes_fn)
		# Load to: state_probs_at_corners_across_all_trees
		load(file=stateprobs_at_corners_fn)
		} # end if runslow==TRUE
	
	
	# Also store parameter estimates
	optim_results_table = NULL
	get_optim_results = TRUE
	if (get_optim_results == TRUE)
		{
		cat("\n\nGetting ML parameter estimates for trees #", start_treenum, "-", end_treenum, ":\n", sep="")
		for (i in start_treenum:end_treenum)
		#for (i in 1:length(newick_fns))
		#for (i in 1:84)
			{
			cat(i, " ", sep="")

			# Get the BioGeoBEARS results for this tree
			# Loads to "res"
			resfn = resfns[i]
			load(resfn)
	
			# Store the parameter estimates
			optim_results_table = rbind(optim_results_table, res$optim_result)
			} # end for (i in 1:length(newick_fns))
		} # end if (get_optim_results == TRUE)

	optim_results_mean = colMeans(optim_results_table)
	optim_results_sd = apply(X=optim_results_table, MARGIN=2, FUN=sd)
	
	# These results are atomic vectors, convert to data.frames as in optimx
	tmp_colnames = names(optim_results_mean)
	optim_results_mean = data.frame(matrix(data=optim_results_mean, nrow=1))
	names(optim_results_mean) = tmp_colnames
	optim_results_sd = data.frame(matrix(data=optim_results_sd, nrow=1))
	names(optim_results_sd) = tmp_colnames
	
	# Return results
	results_on_multiple_trees = list()
	
	
	# Add the state probabilities, if you get those...
	if ((run_or_collect == "collect") || (run_or_collect == "both"))
		{
		results_on_multiple_trees$state_probs_at_nodes_across_all_trees = state_probs_at_nodes_across_all_trees
		results_on_multiple_trees$state_probs_at_corners_across_all_trees = state_probs_at_corners_across_all_trees
		}

	# Filenames
	results_on_multiple_trees$newick_fns = newick_fns
	results_on_multiple_trees$geog_fns = geog_fns
	results_on_multiple_trees$resfns = resfns
	
	# Optim results (fast)
	results_on_multiple_trees$optim_results_table = optim_results_table
	results_on_multiple_trees$optim_results_mean = optim_results_mean
	results_on_multiple_trees$optim_results_sd = optim_results_sd
	
	
	if (plot_params == TRUE)
		{
		plot_params_from_multiple_trees(optim_results_table, BioGeoBEARS_run_object, optimx2012=FALSE)
		}
	
	return(results_on_multiple_trees)
	}



plot_params_from_multiple_trees <- function(optim_results_table, BioGeoBEARS_run_object, optimx2012=FALSE)
	{
	if (optimx2012 == TRUE)
		{
		stop("ERROR in plot_params_from_multiple_trees(): This function is not yet implemented for optimx 2012")
		}
	
	# Find the colnums and names of the free parameters
	value_TF = names(optim_results_mean) == "value"
	tmp_colnums = 1:length(value_TF)
	LnL_colnum = tmp_colnums[value_TF]
	nparams = LnL_colnum - 1

	# parameter names
	TF = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table$type == "free"
	param_names = rownames(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table)[TF]

	numgraphs = nparams+1
	#par(mfrow=c(numgraphs, 1))

	# Get min/max
	minval = 0
	maxval = max(optim_results_table[,1:nparams], na.rm=TRUE)
	maxx = max(pretty(c(minval, maxval)))
	for (j in 1:nparams)
		{
		titletxt = paste("distribution of ", param_names[j], " over ", length(newick_fns), " trees", sep="")
		hist(optim_results_table[,j], main=titletxt, xlim=c(minval, maxx), xlab=param_names[j], breaks=30)
		}
	titletxt = paste("distribution of LnL over ", length(newick_fns), " trees", sep="")
	hist(optim_results_table[,LnL_colnum], main=titletxt, xlab="LnL", breaks=30)
	
	# Scatterplots
	for (i in 1:(nparams-1))
		{
		for (j in 2:nparams)
			{
			plot(x=optim_results_table[,i], y=optim_results_table[,j], xlim=c(minval, maxx), ylim=c(minval, maxx), xlab=param_names[i], ylab=param_names[j], main="")
			}
		}
	}


# Summarize the range/state probabilities on a master tree
summarize_stateprobs_on_master_tree <- function(master_tree_fn, state_probs_at_nodes_across_all_trees, state_probs_at_corners_across_all_trees, plotflag=FALSE)
	{
	defaults='
	master_tree_fn = "Finalgent_last500_MCC_378tips.newick"
	'
	
	# Number of states is the number of columns in state probabilities,
	# minus the OTU names in the last column
	numstates = ncol(state_probs_at_nodes_across_all_trees) - 1

	#######################################################
	# Take the state probabilities from the sample of trees
	# that were run on BioGeoBEARS, and put them on the 
	# nodes of your master / (bifurcating) consensus tree
	#######################################################

	# Load the master_tree
	master_tree = read.tree(file=master_tree_fn)
	
	if (plotflag)
		{
		plot(master_tree, show.tip.label=FALSE)
		axisPhylo()
		title("Master tree")
		length(master_tree$tip.label)
		}




	# Now, here's the tricky bit.  
	# The master tree has 378 living tips.
	# The sample of Newick trees has 408 tips (378 living, the rest fossils)
	# *IF* the master tree has no extra tips not in the sampled trees, you can do this:


	#######################################################
	# Get a list of all OTUs
	#######################################################
	# Check that master_tree has nothing outside of this list
	OTUs_list_of_lists = unique(state_probs_at_nodes_across_all_trees[,ncol(state_probs_at_nodes_across_all_trees)])
	length(OTUs_list_of_lists)
	# Paste them all together
	OTUs_non_unique = paste(OTUs_list_of_lists, collapse=",", sep="")
	OTUs_unique = sort(unique(strsplit(OTUs_non_unique, split=",")[[1]]))
	length(OTUs_unique)

	TF = master_tree$tip.label %in% OTUs_unique
	if (sum(TF) < length(TF))
		{
		errortxt = "ERROR: Your master_tree has tips that are NOT in your sampled tree tips.\nYou should reduce the master_tree with drop.tip.\nTips to drop are:\n\n"
		cat(errortxt)
		cat(master_tree$tip.label[TF], sep="\n")
		cat("\n\n")
		}

	# Basically, we want to find the SMALLEST clade including everything of interest
	# clade_size = # of commas + 1
	library(stringr)	# for str_count
	numcommas_all_OTUs = str_count(string=OTUs_list_of_lists, pattern=",")
	numOTUs_all_OTUs = numcommas_all_OTUs + 1

	# Translate the full list of names, to the list of names without fossils
	OTUs_not_in_master_tree_TF = (OTUs_unique %in% master_tree$tip.label == FALSE)
	sum(OTUs_not_in_master_tree_TF)
	OTUs_not_in_master_tree = OTUs_unique[OTUs_not_in_master_tree_TF]

	# Make a translated list of unique OTU lists
	# Remove all OTUs not in the master tree
	OTUs_list_of_lists_reduced = OTUs_list_of_lists
	for (i in 1:length(OTUs_not_in_master_tree))
		{
		OTUs_list_of_lists_reduced = gsub(pattern=OTUs_not_in_master_tree[i], replacement="", x=OTUs_list_of_lists_reduced)
		}
	
	# Remove all resulting double-commas, beginning/ending commas
	OTUs_list_of_lists_reduced = remove_multiple_commas_from_strings_in_list(list_of_strings=OTUs_list_of_lists_reduced)
	

	# Count the remaining commas
	numcommas_reduced_OTUs = str_count(string=OTUs_list_of_lists_reduced, pattern=",")
	numOTUs_reduced_OTUs = numcommas_reduced_OTUs + 1

	# Check that the OTU numbers match in your nodes list vs
	# your master_tree
	if (max(numOTUs_reduced_OTUs) != length(master_tree$tip.label) )
		{
		errortxt = "\n\nERROR: max(numOTUs_reduced_OTUs) != length(master_tree$tip.label)\nThese should match!\n\n"
		cat(errortxt)
		cat("\nmax(numOTUs_reduced_OTUs)=", max(numOTUs_reduced_OTUs), sep="")
		cat("\nlength(master_tree$tip.label)=", length(master_tree$tip.label), sep="")
		stop("\nStopping on error. Check your inputs.\n")
		}



	# OK, now we have a list of reduced OTU lists corresponding to each OTUs_list_of_lists
	# Basically, we take the smallest match, identify the rows, and add them up
	length(OTUs_list_of_lists)
	length(OTUs_list_of_lists_reduced)

	# Get the OTU lists for the master_tree
	master_trtable = prt(master_tree, printflag=FALSE, get_tipnames=TRUE)
	master_tr_OTUs = master_trtable$tipnames

	# Go through the unique clades in the MCMC tree sample
	numnodes_master_tr = length(master_tr_OTUs)
	nummatches_at_node = rep(NA, times=numnodes_master_tr)
	stateprobs_nodes_master_tr = matrix(NA, nrow=numnodes_master_tr, ncol=numstates)
	stateprobs_corners_master_tr = matrix(NA, nrow=numnodes_master_tr, ncol=numstates)


	cat("\n\nProcessing master tree node #", sep="")
	for (i in 1:length(master_tr_OTUs))
	#for (i in 1:1)
		{
		cat(i, ",", sep="")
	
		master_tr_OTU = master_tr_OTUs[i]
		master_matches_in_OTUs_reduced_list_TF = OTUs_list_of_lists_reduced == master_tr_OTU

		# Better not be 0!!!
		if (sum(master_matches_in_OTUs_reduced_list_TF) == 0)
			{
			cat("\n\nERROR on node #", i, ": 0 matchs of clade OTUs  to OTUs_list_of_lists_reduced!!\n\n", sep="")
		
			cat("\nClade OTUs are:\n\n")
			cat(master_tr_OTU)
			cat("\n")
			#print(minTF)
			#print(sum(minTF))
			#print(OTUs_list_of_lists[master_matches_in_OTUs_reduced_list_TF][minTF])
			#stop("\n\nStopped on error.")
		
			cat("\n\nERROR on node #", i, ": 0 matchs of clade OTUs  to OTUs_list_of_lists_reduced!!", sep="")
			cat("\n\nThis may not be unusual, if you have some clades in your master tree which are poorly-supported\n(usually your master_tree is a consensus tree (should be strictly bifurcating) or\nMCC tree (which is just one of the trees from the BEAST MCMC sample, and may well have\nsome poorly supported clades), and you are only looking at ", length(newick_fns), " sampled trees\nfrom the MCMC sample.\n\n", sep="")
			cat("\n\nWe'll just set the state probabilities of this node to NA.\n\n", sep="")
		
			nummatches_at_node[i] = 0
			stateprobs_nodes_master_tr[i, ] = NA
			stateprobs_corners_master_tr[i, ] = NA
			next()
			}


		# Get the smallest one
		minval = min(numOTUs_all_OTUs[master_matches_in_OTUs_reduced_list_TF])
		minTF = numOTUs_all_OTUs[master_matches_in_OTUs_reduced_list_TF] == minval



	
		# Better be 1!!!
		if (sum(minTF) != 1)
			{
			cat("\n\nERROR on node #", i, ": More than one match of clade OTUs '", master_tr_OTU, "' to OTUs_list_of_lists_reduced!!\n\n", sep="")
			print(minTF)
			print(sum(minTF))
			print(OTUs_list_of_lists[master_matches_in_OTUs_reduced_list_TF][minTF])
			stop("\n\nStopped on error:  More than one match of clade OTUs to OTUs_list_of_lists_reduced")
			}
	
		# Identify the rows for matching nodes in the MCMC tree samples
		OTUs_of_a_node_in_the_fossil_trees = OTUs_list_of_lists[master_matches_in_OTUs_reduced_list_TF][minTF]
	
		# Get the probabilities for the nodes
		TF = state_probs_at_nodes_across_all_trees[,ncol(state_probs_at_nodes_across_all_trees)] == OTUs_of_a_node_in_the_fossil_trees
		temp_probsmat = matrix(state_probs_at_nodes_across_all_trees[TF, 1:numstates], nrow=sum(TF))
		state_probs_for_matching_node = matrix(apply(X=temp_probsmat, MARGIN=2, FUN=as.numeric), nrow=sum(TF))
	
		nummatches_at_node[i] = nrow(state_probs_for_matching_node)
		colSums_state_probs_for_matching_node = colSums(state_probs_for_matching_node, na.rm=TRUE)
	
		stateprobs_nodes_master_tr[i, ] = colSums_state_probs_for_matching_node / sum(colSums_state_probs_for_matching_node)


		# Get the probabilities for the corners
		TF = state_probs_at_corners_across_all_trees[,ncol(state_probs_at_corners_across_all_trees)] == OTUs_of_a_node_in_the_fossil_trees
		temp_probsmat = matrix(state_probs_at_corners_across_all_trees[TF, 1:numstates], nrow=sum(TF))
		state_probs_for_matching_corner = matrix(apply(X=temp_probsmat, MARGIN=2, FUN=as.numeric), nrow=sum(TF))
	
		nummatches_at_node[i] = nrow(state_probs_for_matching_corner)
		colSums_state_probs_for_matching_corner = colSums(state_probs_for_matching_corner, na.rm=TRUE)
	
		stateprobs_corners_master_tr[i, ] = colSums_state_probs_for_matching_corner / sum(colSums_state_probs_for_matching_corner)

		} # end for (i in 1:length(master_tr_OTUs))

	head(stateprobs_nodes_master_tr)
	head(stateprobs_corners_master_tr)

	rowSums(stateprobs_nodes_master_tr)
	rowSums(stateprobs_corners_master_tr)
	
	stateprobs_list = list()
	stateprobs_list$stateprobs_nodes_master_tr = stateprobs_nodes_master_tr
	stateprobs_list$stateprobs_corners_master_tr = stateprobs_corners_master_tr
	stateprobs_list$nummatches_at_node = nummatches_at_node
	
	return(stateprobs_list)
	}


remove_multiple_commas_from_strings_in_list <- function(list_of_strings)
	{
	# Convert variable name
	OTUs_list_of_lists_reduced = list_of_strings
	
	# Remove all resulting extraneous commas
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,,", replacement=",", x=OTUs_list_of_lists_reduced)
	OTUs_list_of_lists_reduced = gsub(pattern=",,", replacement=",", x=OTUs_list_of_lists_reduced)

	# Remove commas at the end of strings
	lengths = sapply(X=OTUs_list_of_lists_reduced, FUN=nchar)
	names(lengths) = NULL
	last_chars = mapply(FUN=substr, OTUs_list_of_lists_reduced, start=lengths, stop=lengths)
	names(last_chars) = NULL
	ending_comma_TF = last_chars == ","

	stops = -1 + lengths[ending_comma_TF]
	starts = rep(1, times=length(stops))

	no_ending_commas = mapply(FUN=substr, OTUs_list_of_lists_reduced[ending_comma_TF], start=starts, stop=stops)
	names(no_ending_commas) = NULL
	no_ending_commas
	OTUs_list_of_lists_reduced[ending_comma_TF] = no_ending_commas


	# Remove commas at the beginning of strings
	first_chars = sapply(X=OTUs_list_of_lists_reduced, FUN=substr, start=1, stop=1)
	names(first_chars) = NULL
	starting_comma_TF = first_chars == ","
	sum(starting_comma_TF)
	#cbind(first_chars, starting_comma_TF)[500:750,]

	starts = rep(2, times=sum(starting_comma_TF))
	stops = sapply(X=OTUs_list_of_lists_reduced[starting_comma_TF], FUN=nchar)

	no_starting_commas = mapply(FUN=substr, OTUs_list_of_lists_reduced[starting_comma_TF], start=starts, stop=stops)
	names(no_starting_commas) = NULL
	OTUs_list_of_lists_reduced[starting_comma_TF] = no_starting_commas

	#OTUs_list_of_lists_reduced[500:750]

	sum(str_count(string=OTUs_list_of_lists_reduced, pattern="fossil"))
	
	return(OTUs_list_of_lists_reduced)
	}


make_BioGeoBEARS_manytrees_results_object <- function(BioGeoBEARS_run_object, master_tree_fn, geogfn, stateprobs_list, optim_results_mean)
	{
	BioGeoBEARS_results_object_manytrees = BioGeoBEARS_run_object
	BioGeoBEARS_results_object_manytrees$trfn = master_tree_fn
	BioGeoBEARS_results_object_manytrees$geog = geogfn
	
	bears_output_manytrees = list()
	bears_output_manytrees$inputs = BioGeoBEARS_results_object_manytrees
	bears_output_manytrees$spPmat_inputs = NULL
	
	# Store the input BioGeoBEARS_model_object; then modify
	bears_output_manytrees$outputs = BioGeoBEARS_run_object$BioGeoBEARS_model_object
	params = BioGeoBEARS_model_object_to_est_params(BioGeoBEARS_run_object$BioGeoBEARS_model_object)
	
	bears_output_manytrees$optim_result = optim_results_mean
	
	bears_output_manytrees$ML_marginal_prob_each_state_at_branch_top_AT_node = stateprobs_list$stateprobs_nodes_master_tr
	bears_output_manytrees$ML_marginal_prob_each_state_at_branch_bottom_below_node = stateprobs_list$stateprobs_corners_master_tr
	bears_output_manytrees$total_loglikelihood = get_LnL_from_optim_result(optim_results_mean, use_optimx=BioGeoBEARS_run_object$use_optimx)
	
	return(bears_output_manytrees)
	}


