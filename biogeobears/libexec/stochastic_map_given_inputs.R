#######################################################
# Stochastic mapping on a non-stratified analysis
#######################################################
# 
# This function does stochastic mapping on non-stratified analysis, or 
# calls stochastic_mapping_on_stratified() if the input is stratified.
#  
# stochastic_mapping_inputs = get_inputs_for_stochastic_mapping_from_results_object(res=resDEC)
# stochastic_mapping_results = stochastic_map_given_inputs(stochastic_mapping_inputs)
# 
# For stochastic mapping on stratified analysis:
# 
# stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping_stratified(res=resDEC)
# stochastic_mapping_results = stochastic_mapping_on_stratified(res=resDEC, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list)
# 


stochastic_map_given_inputs <- function(stochastic_mapping_inputs, piecenum=NULL, maxtries=40000, seedval=as.numeric(Sys.time()), include_null_range, master_nodenum_toPrint=0)
	{
	#######################################################
	# running_stochastic_mapping, given a results object
	#######################################################
	defaults='
	res
	rootedge=FALSE
	statenum_bottom_root_branch_1based=NULL
	printlevel=1
	stratified=FALSE
	min_branchlength=0.000001
	
	stochastic_mapping_inputs=stochastic_mapping_inputs_list
	piecenum=NULL
	maxtries=40000
	seedval=12346	
	include_null_range=TRUE
	'
	
	# Set the seed
	set.seed(seedval)
	
	# Necessary setting to avoid getting numbers etc. in the stochastic mapping output
	options(stringsAsFactors=FALSE)
		
	# Load inputs 
	rootedge = stochastic_mapping_inputs$rootedge
	statenum_bottom_root_branch_1based = stochastic_mapping_inputs$statenum_bottom_root_branch_1based
	printlevel = stochastic_mapping_inputs$printlevel
	stratified = stochastic_mapping_inputs$stratified
	
	if (is.null(stratified) == TRUE)
		{
		stop("ERROR in stochastic_map_given_inputs(); stochastic_mapping_inputs$stratified is NULL")
		}
	
	
	min_branchlength = stochastic_mapping_inputs$min_branchlength
	cluster_already_open = stochastic_mapping_inputs$cluster_already_open
	res = stochastic_mapping_inputs$res
	Qmat = stochastic_mapping_inputs$Qmat
	COO_weights_columnar = stochastic_mapping_inputs$COO_weights_columnar
	Rsp_rowsums = stochastic_mapping_inputs$Rsp_rowsums

	tipranges = stochastic_mapping_inputs$tipranges
	areas = stochastic_mapping_inputs$areas
	state_indices_0based = stochastic_mapping_inputs$state_indices_0based
	states_list = state_indices_0based
	ranges_list = stochastic_mapping_inputs$ranges_list
	numstates = stochastic_mapping_inputs$numstates

	if (include_null_range == TRUE)
		{
		numstates_during_cladogenesis = numstates - 1
		} else {
		numstates_during_cladogenesis = numstates - 0
		}
	
	
	
	if (stratified == FALSE)
		{
		trtable = stochastic_mapping_inputs$trtable
		tr = stochastic_mapping_inputs$tr
		phy2 = stochastic_mapping_inputs$phy2
		
		independent_likelihoods_on_each_branch = stochastic_mapping_inputs$independent_likelihoods_on_each_branch

		
		} else {
		subtable_rowsTF = stochastic_mapping_inputs$master_table_timeperiod_i$piecenum == piecenum
		subtable_rownums = (1:nrow(stochastic_mapping_inputs$master_table_timeperiod_i))[subtable_rowsTF]
		
 		# Convert master_table_timeperiod_i to trtable
		trtable = stochastic_mapping_inputs$master_table_timeperiod_i[subtable_rownums,]
		# trtable parts to use/change are ONLY the parts corresponding to this particular treepiece
		#trtable_rows_correct_pieces_TF = trtable_all_pieces$piecenum == piecenum
		#trtable = trtable_all_pieces[trtable_rows_correct_pieces_TF,]
		
		independent_likelihoods_on_each_branch = stochastic_mapping_inputs$independent_likelihoods_by_tree_piece_for_timeperiod_i[[piecenum]]
		
		phy = stochastic_mapping_inputs$tree_sections_list$return_pieces_list[[piecenum]]
		phy2 = reorder(phy, "pruningwise") # Do this, 
		} # END if (stratified == FALSE)


	# Basic tree info
	ntips = length(phy2$tip.label)
	num_internal_nodes = phy2$Nnode
	tipnums = 1:ntips
	root_nodenum = ntips+1
	nodenums = root_nodenum:(ntips+num_internal_nodes)
	
	
	


	if (stratified == FALSE)
		{
		# Add a column for the sampled node states
		sampled_states_AT_nodes = rep(NA, nrow(trtable))
		sampled_states_AT_brbots = rep(NA, nrow(trtable))

		# Add the right and left descendant node numbers
		leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(phy2)

		left_desc_nodes = rep(NA, nrow(trtable))
		right_desc_nodes = rep(NA, nrow(trtable))

		# dcorner = descendant corner (i.e. right after speciation)
		samp_LEFT_dcorner = rep(NA, nrow(trtable))
		samp_RIGHT_dcorner = rep(NA, nrow(trtable))

		# Events on branches
		anagenetic_events_txt_below_node = rep(NA, nrow(trtable))

		trtable = cbind(trtable, sampled_states_AT_nodes, sampled_states_AT_brbots, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner, anagenetic_events_txt_below_node)
		
		# This works, the reverse does not -- 2015-04-06
		# (LR for tree plots is opposite for the internal structure)
		trtable$left_desc_nodes[nodenums] = leftright_nodes_matrix$right
		trtable$right_desc_nodes[nodenums] = leftright_nodes_matrix$left
		trtable[nodenums,]
		} else {
		# Stratified analysis: 
		# You've already got all these columns
		junk = 1
		} # END if (stratified == FALSE)


	#returned_mats1 = get_Qmat_COOmat_from_BioGeoBEARS_run_object(default_BioGeoBEARS_run_object)
	#returned_mats1

# 	returned_mats2 = get_Qmat_COOmat_from_res(res)
# 	returned_mats2
# 
# 	# Extract output
# 	Qmat = returned_mats2$Qmat
# 	COO_weights_columnar = returned_mats2$COO_weights_columnar
# 	Rsp_rowsums = returned_mats2$Rsp_rowsums

	# Calculate the likelihood P((left_state,right_state)|anc_state)
	# for each scenario (unconstrained)
	# Note:
	# COO_weights_columnar indices are 0-based, with no null_range
	# So, add 2 to get it so that e.g. state 0 = state 2 = Kauai
	#
	# Or, add 1 to get the 1based state indexes INSIDE COO_weights_columnar
	# 
	# COO_weights_columnar =
	# ancestral index, left index, right index, conditional
	# probability given ancestral states. (assuming likelihood
	# of descendants is 1)
	# Probabilities of each range-inheritance scenario, conditional
	# on ancestral state (without constraints on Left Branch state)
	#like_LeftRight_given_AncState = COO_weights_columnar[[4]] / (Rsp_rowsums[1+COO_weights_columnar[[1]]])
	#like_LeftRight_given_AncState


	# Calculate the total number of range-inheritance scenarios
	# under the model
	# (this is the number of scenarios with weight > 0)
	# (weight per event/ sum(weights) = prob. per event)
	num_scenarios = length(COO_weights_columnar[[1]])



	#######################################################
	#######################################################
	# THIS IS AN UPPASS FROM THE TIPS TO THE ROOT
	#######################################################
	#######################################################

	# Check to make sure you have the necessary inputs
	if (exists("COO_weights_columnar") == FALSE)
		{
		stop("\nERROR_A: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
		}
	if (exists("Rsp_rowsums") == FALSE)
		{
		stop("\nERROR_B: calc_loglike_sp requires 'COO_weights_columnar', 'Rsp_rowsums', and cppSpMethod==3 for marginal ancestral state estimations.\n")
		}

	##########################################################
	# 1. Sample a state at the root of the tree or subtree (if it doesn't exist)
	# 2. Given the root state, calculate uppass probabilities to the corners
	# 3. Multiply the corner uppass probs by the downpass probs
	# 
	##########################################################
	
	
	# This is for sampling just the root node state
	# If there is no rootedge, and if the starting state
	# at the root node is not determined
	
	# If stratified analysis, no root edge, no pre-determined root state
	TF1 = ((stratified == TRUE) && (rootedge == FALSE) && (is.na(trtable$sampled_states_AT_nodes[root_nodenum])) && (trtable$node.type[root_nodenum] == "root"))

	# If non-stratified analysis, whether or not there is a root edge, no pre-determined root state
	#TF2 = ((stratified == FALSE) && (is.na(trtable$sampled_states_AT_nodes[root_nodenum])) && (trtable$node.type[root_nodenum] == "root"))
	TF2 = (stratified == FALSE)
	#print(TF2)
	if ( TF1 || TF2 )
		{
		# Sample a state at the root
		
		# The global root nodenum will be different than the subtree root nodenum
		if (stratified == TRUE)
			{
			global_rootTF = trtable$node.type == "root"
			# Error check
			if (sum(global_rootTF) != 1)
				{
				stop("\n\nStop ERROR: no global root node in subtree.\n\n")
				}
			global_root_nodenum = trtable$node[global_rootTF]
			
			} else {
			global_root_nodenum = root_nodenum
			} # END if (stratified == TRUE)
		
		# Get the stateprobs at the global root, and sample from them
		probs_branch_top = res$ML_marginal_prob_each_state_at_branch_top_AT_node[global_root_nodenum,]
		statenums = 1:numstates
		statenum_1based = sample(x=statenums, size=1, replace=TRUE, prob=probs_branch_top)
		statenum_1based
		#print(global_root_nodenum)
		#print(round(node_stateprobs,3))
		#print(statenum_1based)
		
		# Store the sampled state at the root
		trtable$sampled_states_AT_nodes[root_nodenum] = statenum_1based
		}

	
	# If stratified analysis, with root edge, no pre-determined root state
	# Simulate up the root branch, if that exists
	# (and if its a stratified analysis)
	#print(stratified)
	if ( (stratified == TRUE) && (rootedge == TRUE))
		{
		# Find the downpass conditional likelihoods (normalized) that have
		# been pre-calculated
		#res$condlikes
		#res$inputs$master_table
		subtable_rootTF = trtable$SUBnode.type == "root"
		subtable_rownum = (1:nrow(trtable))[subtable_rootTF]
		
		TF1 = res$inputs$master_table$node == trtable$node[subtable_rownum]
		TF2 = res$inputs$master_table$stratum == trtable$stratum[subtable_rownum]
		TF3 = res$inputs$master_table$piecenum == trtable$piecenum[subtable_rownum]
		TF4 = res$inputs$master_table$piececlass == "subtree"
		TF5 = res$inputs$master_table$SUBnode.type == "root"
		TF = (TF1 + TF2 + TF3 + TF4 + TF5) == 5
		rownums = 1:nrow(res$inputs$master_table)
		rownum = rownums[TF]
		rownum
		downpass_condlikes_at_branch_top = res$condlikes[rownum,]
		downpass_relprobs_at_branch_top = downpass_condlikes_at_branch_top / sum(downpass_condlikes_at_branch_top)
		
		# Now you just need to exponentiate up, given the previous-done 
		# independent likelihoods
		starting_state_1based = trtable$sampled_states_AT_brbots[root_nodenum]
		condprobs_branch_top = rep(0, times=numstates)
		condprobs_branch_bot = rep(0, times=numstates)
		condprobs_branch_bot[starting_state_1based] = 1
		
		# Exponentiate up (well, sorta, exponential pre-calculated)
		#condprobs_branch_top = condprobs_branch_bot %*% independent_likelihoods_by_tree_piece_for_timeperiod_i[[piecenum]]
		branch_length = trtable$SUBedge.length[root_nodenum]
		independent_likelihoods_on_root_branch_of_subtree = expokit_dgpadm_Qmat2(times=branch_length, Qmat=Qmat, transpose_needed=TRUE)
		
		condprobs_branch_top = condprobs_branch_bot %*% independent_likelihoods_on_root_branch_of_subtree
		
		if (include_null_range == TRUE)
			{
			condprobs_branch_top[1] = 0	# zero out the NULL range, since it is impossible in a survivor
			}
		
		# State probabilities at the top of the branch
		probs_branch_top = condprobs_branch_top * downpass_relprobs_at_branch_top
		probs_branch_top = probs_branch_top / sum(probs_branch_top)
		
		

		master_nodenum = res$inputs$master_table$node[rownum]
		if (master_nodenum == master_nodenum_toPrint)
			{
			print("stochastic_map_given_inputs():")
			print("rownum:")
			print(rownum)
			print("master_nodenum_toPrint:")
			print(master_nodenum_toPrint)
			print("condprobs_branch_bot:")
			print(round(condprobs_branch_bot, 3))
			print("downpass_relprobs_at_branch_top:")
			print(round(downpass_relprobs_at_branch_top, 3))
			print("condprobs_branch_top:")
			print(round(condprobs_branch_top, 3))
			print("probs_branch_top:")
			print(round(probs_branch_top, 3))
			}
		
		
		
		
		# Sample the state
		sampled_state_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_branch_top)
		sampled_state_branch_top_1based
		
		# Store the state
		trtable$sampled_states_AT_nodes[subtable_rownum] = sampled_state_branch_top_1based


		#######################################################
		# Stochastic mapping, once the states at branch bottoms
		# and branch tops have been sampled
		# Specifically for root branches of sub-trees!!
		# (left this out the first time -- 2014-05-28_NJM)
		#######################################################
		
		# Stochastic mapping of events on the subtree root branch
		events_table_for_branch_below_subtree_root_node = stochastic_map_branch(nodenum_at_top_of_branch=subtable_rownum, trtable=trtable, Qmat=Qmat, state_indices_0based=state_indices_0based, ranges_list=ranges_list, areas=areas, stratified=stratified, maxtries=maxtries)

		# Store the text representation
		# (extract to table with events_txt_into_events_table() )
		subtree_root_branch_events_txt = events_table_into_txt(events_table_for_branch_below_subtree_root_node)
	
		trtable$anagenetic_events_txt_below_node[subtable_rownum] = subtree_root_branch_events_txt
		# End stochastic mapping on the branches below subtree root
		} # END if ( (stratified == TRUE) && (rootedge == TRUE))
	

	
	#######################################################
	# UPPASS THROUGH THE NODE FROM THE ROOT NODE
	# START FROM A NODE STATE, THEN SIMULATE THE TWO NODE STATES ABOVE,
	# THEN SIMULATE THE BRANCH EVENTS
	#######################################################
	
	# Visit edges in reverse order from the downpass
	edges_to_visit_uppass = seq(from=(num_internal_nodes*2), by=-2, length.out=num_internal_nodes)
	# Since we are going backwards
	#print(edges_to_visit_uppass)
	#print(i)
	#print(j)
	#cat("\n")

	#for (i in edges_to_visit_uppass)
	#j=edges_to_visit_uppass[1]
	for (j in edges_to_visit_uppass)		# Since we are going backwards
		{
		# First edge visited is i
		#print(i)
	
		# Its sister is j 
		#j <- i - 1
		i <- j - 1		# Since we are going backwards

		# Get the node numbers at the tips of these two edges		
		left_desc_nodenum <- phy2$edge[i, 2]
		right_desc_nodenum <- phy2$edge[j, 2]

		# And for the ancestor edge (i or j shouldn't matter, should produce the same result!!!)
		anc <- phy2$edge[i, 1]
		# Store the node number (starting with the root)
		nodenum = anc
	
		# get the correct edges
		left_edge_TF = phy2$edge[,2] == left_desc_nodenum
		right_edge_TF = phy2$edge[,2] == right_desc_nodenum
	
		# Check the branchlength of each edge
		# It's a hook if either branch is super-short
		is_leftbranch_hook_TF = phy2$edge.length[left_edge_TF] < min_branchlength
		is_rightbranch_hook_TF = phy2$edge.length[right_edge_TF] < min_branchlength
		hooknode_TF = (is_leftbranch_hook_TF + is_rightbranch_hook_TF) > 0



		# Get the state at the current node (anc)
		statenum_1based = trtable$sampled_states_AT_nodes[nodenum] 

		#######################################################
		# STOCHASTIC MAPPING OF CLADOGENETIC PROCESS
		#######################################################
		# If it's a hooknode, then copy the node state up
		# (no sampling needed)
		if (hooknode_TF == TRUE)
			{
			# Just copy the node state to the corners
			sampled_split_descendants = list()
			sampled_split_descendants$left_decstate_1based = statenum_1based
			sampled_split_descendants$right_decstate_1based = statenum_1based
			# These will get copied to the table outside of this loop
			} else {
			# If NOT a hooknode (typical), sample a pair of descendant states	
			# Calculate the probability of each range inheritance scenario, 
			# given the chosen root state
			
			# -1 regardless of whether there is a null range (1 to 0 conversion)
			index_Qmat_0based_of_starting_state = statenum_1based - 1
	
			#RCOO_probs_list_given_ancestor = given_a_starting_state_get_prob_of_each_split_scenario(index_Qmat_0based_of_starting_state, COO_weights_columnar, numstates=numstates_during_cladogenesis, include_null_range=TRUE)
	
			#uppass_probs_of_scenarios_given_root_state = RCOO_probs_list_given_ancestor
			#uppass_probs_of_scenarios_given_root_state

			#cbind(COO_weights_columnar[[1]], COO_weights_columnar[[2]], COO_weights_columnar[[3]], uppass_probs_of_scenarios_given_root_state)

			# The downpass probs of each state at each branch
			if (stratified == FALSE)
				{
				# stratified == FALSE
				# *************** perhaps use 
				# left_desc_nodenum and right_desc_nodenum
				left_branch_decnode = trtable$left_desc_nodes[nodenum]
				right_branch_decnode = trtable$right_desc_nodes[nodenum]

				# *************** perhaps use 
				# left_desc_nodenum and right_desc_nodenum
				# NOPE, THIS IS FINE -- 2014-12-29_NJM
				if ( (left_desc_nodenum != left_branch_decnode) | (right_desc_nodenum != right_branch_decnode) )
					{
					print("left_desc_nodenum, right_desc_nodenum")
					cat(left_desc_nodenum, right_desc_nodenum)
					
					print("left_branch_decnode, right_branch_decnode")
					cat(left_branch_decnode, right_branch_decnode)
					
					stop()
					}


				left_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_branch_decnode, ]
				right_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_branch_decnode, ]
				} else {
				# stratified==TRUE
				daughter_nodenums_global = trtable$daughter_nds[nodenum][[1]]
				# names(leftright_nodes_matrix) = c("right", "left")
				left_branch_decnode_global = daughter_nodenums_global[1]
				right_branch_decnode_global = daughter_nodenums_global[2]
				left_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[left_branch_decnode_global, ]
				right_branch_downpass_likes = res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS[right_branch_decnode_global, ]
				
				# Get the LOCAL left/right nodes, for the 
				# sub-table
				left_branch_decnode_TF = trtable$node == left_branch_decnode_global
				left_branch_decnode = (1:nrow(trtable))[left_branch_decnode_TF]

				right_branch_decnode_TF = trtable$node == right_branch_decnode_global
				right_branch_decnode = (1:nrow(trtable))[right_branch_decnode_TF]
				
				left_branch_decnode
				right_branch_decnode
				} # END if (stratified == FALSE)
				
			#sampled_split_descendants = sample_split_scenario(COO_weights_columnar, uppass_probs_of_scenarios_given_root_state, left_branch_downpass_likes, right_branch_downpass_likes)
			
			# Input the ancestral state into the state probabilities at the current node
			statenum_1based
			probs_ancstate = rep(0, length(left_branch_downpass_likes))
			probs_ancstate[statenum_1based] = 1
			probs_ancstate
			
			# NJM -- for bug checking...
			if (nodenum == 21)
				{
				printflag = TRUE
				} else {
				printflag = FALSE
				}
			
			# OLD (2014-05-ish)
			#sampled_split_descendants = sample_split_scenario2(COO_weights_columnar, probs_ancstate, left_branch_downpass_likes, right_branch_downpass_likes, sample_which="both", return_prob_each_split_scenario=FALSE, include_null_range=TRUE, Rsp_rowsums=NULL, numstates_wo_null=NULL, printflag=printflag)
			
			# NEW (2015-01-10)
			sample_uppass_res = sample_uppass_split_scenario_given_probs_ancstate(probs_ancstate=probs_ancstate, COO_weights_columnar=COO_weights_columnar, numstates=numstates, include_null_range=include_null_range, left_branch_downpass_likes=left_branch_downpass_likes, right_branch_downpass_likes=right_branch_downpass_likes, Rsp_rowsums=NULL)

			} # END if (hooknode_TF == TRUE)
		# OLD 
		#left_decstate_1based = sampled_split_descendants$left_decstate_1based
		#right_decstate_1based = sampled_split_descendants$right_decstate_1based
		# NEW 
		left_decstate_1based = as.numeric(sample_uppass_res$left_decstate_1based)
		right_decstate_1based = as.numeric(sample_uppass_res$right_decstate_1based)
		
		#print(c("oldLeft", "oldRight", "newLeft", "newRight"))
		#print(c(sampled_split_descendants$left_decstate_1based, sampled_split_descendants$right_decstate_1based, sample_uppass_res$left_decstate_1based, sample_uppass_res$right_decstate_1based))
		
		# Put these into the trtable
		trtable$samp_LEFT_dcorner[nodenum] = left_decstate_1based
		trtable$samp_RIGHT_dcorner[nodenum] = right_decstate_1based
	
		# And put them in as the sampled states at branch bottoms for the appropriate
		# descendant nods
		# left_branch_decnode and right_branch_decnode are LOCAL
		# for the SUBTREE
		trtable$sampled_states_AT_brbots[left_branch_decnode] = left_decstate_1based
		trtable$sampled_states_AT_brbots[right_branch_decnode] = right_decstate_1based
	
	
		#######################################################
		# STOCHASTIC MAPPING OF ANAGENETIC PROCESS
		#######################################################
		# Now, evolution ALONG branches
		#independent_likelihoods_on_each_branch = calc_independent_likelihoods_on_each_branch(phy2, Qmat, cluster_already_open=NULL, Qmat_is_sparse=FALSE)
		# Steps:
		# a. Given a state at a corner, calculate the conditional probabilities
		#    of states at the branch top.
		# b. Multiply these by the saved downpass probabilities
		# c. Sample from this distribution, & store at the nodes at the top
	
		# Initialize the starting probabilities at branch bottoms
		# (setting the P(known sampled state) to equal 1!!)
		condprobs_Left_branch_top = rep(0, times=numstates)
		condprobs_Right_branch_top = rep(0, times=numstates)

		condprobs_Left_branch_bot = rep(0, times=numstates)
		condprobs_Right_branch_bot = rep(0, times=numstates)
		condprobs_Left_branch_bot[left_decstate_1based] = 1
		condprobs_Right_branch_bot[right_decstate_1based] = 1
	
		# Dense matrix exponentiation, which has been done already!
		TF2 = ( (length(cluster_already_open)==1) && (cluster_already_open==FALSE) )
		if (is.null(cluster_already_open) || (TF2))
			{
			# Relative probabilities of states at the top of left branch
			condprobs_Left_branch_top = condprobs_Left_branch_bot %*% independent_likelihoods_on_each_branch[,,i]
					
			# Relative probabilities of states at the top of right branch
			condprobs_Right_branch_top = condprobs_Right_branch_bot %*% independent_likelihoods_on_each_branch[,,j]
			} else {
		
			# Here, the independent_likelihoods_on_each_branch are stored in a list of matrices
			# Relative probabilities of states at the top of left branch
			condprobs_Left_branch_top = condprobs_Left_branch_bot %*% independent_likelihoods_on_each_branch[[i]]
					
			# Relative probabilities of states at the top of right branch
			condprobs_Right_branch_top = condprobs_Right_branch_bot %*% independent_likelihoods_on_each_branch[[j]]
			} # END if (is.null(cluster_already_open))

		# zero out the NULL range, since it is impossible in a survivor
		if (include_null_range == TRUE)
			{
			condprobs_Left_branch_top[1] = 0
			condprobs_Right_branch_top[1] = 0
			} # END if (include_null_range == TRUE)


		# Get the probabilities at the branch tops for the two branches above the node
		# under consideration
		# In non-stratified -- these are just the node numbers
		# In stratified analysis:
		# Here, left_desc_nodenum & right_desc_nodenum are for the SUBTREE
		if (stratified == FALSE)
			{
			# OK, now multiply the UPPASS and DOWNPASS probabilities
			probs_Left_branch_top = condprobs_Left_branch_top * res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[left_desc_nodenum,]
			probs_Right_branch_top = condprobs_Right_branch_top * res$relative_probs_of_each_state_at_branch_top_AT_node_DOWNPASS[right_desc_nodenum,]
			
			# In case users want to trace what's going on
			left_desc_nodenum_global = left_desc_nodenum
			left_desc_nodenum_global = right_desc_nodenum
			} else {
			# When stratified == TRUE, we have to dig up the corresponding rows of res
			# to get the correct downpass condlikes
			
			# Left node
			TF1 = res$inputs$master_table$node == trtable$node[left_desc_nodenum]
			TF2 = res$inputs$master_table$stratum == trtable$stratum[left_desc_nodenum]
			TF3 = res$inputs$master_table$piecenum == trtable$piecenum[left_desc_nodenum]
			TF4 = res$inputs$master_table$piececlass == "subtree"
			left_desc_nodenum_global_TF = (TF1 + TF2 + TF3 + TF4) == 4
			left_desc_nodenum_global = (1:nrow(res$condlikes))[left_desc_nodenum_global_TF]
			
			# Right node
			TF1 = res$inputs$master_table$node == trtable$node[right_desc_nodenum]
			TF2 = res$inputs$master_table$stratum == trtable$stratum[right_desc_nodenum]
			TF3 = res$inputs$master_table$piecenum == trtable$piecenum[right_desc_nodenum]
			TF4 = res$inputs$master_table$piececlass == "subtree"
			right_desc_nodenum_global_TF = (TF1 + TF2 + TF3 + TF4) == 4
			right_desc_nodenum_global = (1:nrow(res$condlikes))[right_desc_nodenum_global_TF]

			# OK, now multiply the UPPASS and DOWNPASS probabilities
			probs_Left_branch_top = condprobs_Left_branch_top * res$condlikes[left_desc_nodenum_global,]
			probs_Right_branch_top = condprobs_Right_branch_top * res$condlikes[right_desc_nodenum_global,]
			} # END if (stratified == FALSE)
	
		# Normalize by sum so they add to 1
		probs_Left_branch_top = probs_Left_branch_top / sum(probs_Left_branch_top)
		probs_Right_branch_top = probs_Right_branch_top / sum(probs_Right_branch_top)

# 		print("left_desc_nodenum_global:")
# 		print(left_desc_nodenum_global)
# 		print("right_desc_nodenum_global:")
# 		print(right_desc_nodenum_global)
		
		
		#print("Checking:")
		#print(probs_Left_branch_top)
		#print(probs_Right_branch_top)

		#for (zzz in 1:20)
		#	{
		# Sample states at the top of the two descendant branches
		# (these are conditional on the uppass probs, conditional on the 
		#  present node, AND on the saved downpass probs
		sampled_state_Left_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_Left_branch_top)
		sampled_state_Right_branch_top_1based = sample(x=1:numstates, size=1, replace=TRUE, prob=probs_Right_branch_top)
	
		# Store these states
		trtable$sampled_states_AT_nodes[left_branch_decnode] = sampled_state_Left_branch_top_1based

		trtable$sampled_states_AT_nodes[right_branch_decnode] = sampled_state_Right_branch_top_1based



		if (stratified == FALSE)
			{
			master_nodenum = left_desc_nodenum
			} else {
			master_nodenum = res$inputs$master_table$node[left_desc_nodenum_global]
			} # END if (stratified == FALSE)
		if (master_nodenum == master_nodenum_toPrint)
			{
			cat("\n")
			print("left_branch_decnode:")
			print(left_branch_decnode)
			print("condprobs_Left_branch_bot:")
			print(round(condprobs_Left_branch_bot,3))
			print("condprobs_Left_branch_top:")
			print(round(condprobs_Left_branch_top,3))
			print("res$condlikes[left_desc_nodenum_global,]:")
			print(round(res$condlikes[left_desc_nodenum_global,],3))
			print("probs_Left_branch_top:")
			print(round(probs_Left_branch_top,3))
			print("sampled_state_Left_branch_top_1based:")
			print(round(sampled_state_Left_branch_top_1based,3))
			sampled_stateprobs = rep(0, numstates)
			sampled_stateprobs[sampled_state_Left_branch_top_1based] = 1
			print("sampled_stateprobs:")
			print(sampled_stateprobs)
			cat("\n")
			}


		if (stratified == FALSE)
			{
			master_nodenum = right_desc_nodenum
			} else {
			master_nodenum = res$inputs$master_table$node[right_desc_nodenum_global]
			} # END if (stratified == FALSE)
		if (master_nodenum == master_nodenum_toPrint)
			{
			cat("\n")
			print("right_branch_decnode:")
			print(right_branch_decnode)
			print("condprobs_Right_branch_bot:")
			print(round(condprobs_Right_branch_bot,3))
			print("condprobs_Right_branch_top:")
			print(round(condprobs_Right_branch_top,3))
			print("res$condlikes[right_desc_nodenum_global,]:")
			print(round(res$condlikes[right_desc_nodenum_global,],3))
			print("probs_Right_branch_top:")
			print(round(probs_Right_branch_top,3))
			print("sampled_state_Right_branch_top_1based:")
			print(round(sampled_state_Right_branch_top_1based,3))
			sampled_stateprobs = rep(0, numstates)
			sampled_stateprobs[sampled_state_Right_branch_top_1based] = 1
			print("sampled_stateprobs:")
			print(sampled_stateprobs)
			cat("\n")
			}
		




		
		#txt = paste0("zzz: ", zzz, ", L: ", left_decstate_1based, "->", sampled_state_Left_branch_top_1based, "; R: ", right_decstate_1based, "->", sampled_state_Right_branch_top_1based)
		#cat("\n")
		#cat(txt)
		#} # END zzz
		#cat("\n")
	
		# CHECK LEFT-RIGHT NODE NUMBERING
		check_left_vs_right_numbering = FALSE	
		# NOTE: FOR TREE ITERATIONS, YOU HAVE TO SWITCH
		# LEFT AND RIGHT, WHICH IS DIFFERENT THAN
		# FOR GRAPHING:

		# # Add the right and left descendant node numbers
		# leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(phy2)
		# 
		# left_desc_nodes = rep(NA, nrow(trtable))
		# right_desc_nodes = rep(NA, nrow(trtable))
		# 
		# # dcorner = descendant corner (i.e. right after speciation)
		# samp_LEFT_dcorner = rep(NA, nrow(trtable))
		# samp_RIGHT_dcorner = rep(NA, nrow(trtable))
		# 
		# trtable = cbind(trtable, left_desc_nodes, right_desc_nodes, samp_LEFT_dcorner, samp_RIGHT_dcorner)
		# trtable$left_desc_nodes[nodenums] = leftright_nodes_matrix$right
		# trtable$right_desc_nodes[nodenums] = leftright_nodes_matrix$left
		# trtable[nodenums,]

		if (check_left_vs_right_numbering == TRUE)
			{
			if (left_branch_decnode <= ntips)
				{
				print("Left tip:")
				print(left_branch_decnode)
				print(left_desc_nodenum)
				print(sampled_state_Left_branch_top_1based)
				print(probs_Left_branch_top)
				} # END if (left_branch_decnode <= ntips)

			if (right_branch_decnode <= ntips)
				{
				print("Right tip:")
				print(right_branch_decnode)
				print(right_desc_nodenum)
				print(sampled_state_Right_branch_top_1based)
				print(probs_Right_branch_top)
				} # END if (right_branch_decnode <= ntips)
			} # END if (check_left_vs_right_numbering == TRUE)
	
	
		#######################################################
		# Stochastic mapping, once the states at branch bottoms
		# and branch tops have been sampled
		#######################################################
		
		# Stochastic mapping of events on the left branch
		events_table_for_branch_below_Left_node = stochastic_map_branch(nodenum_at_top_of_branch=left_branch_decnode, trtable=trtable, Qmat=Qmat, state_indices_0based=state_indices_0based, ranges_list=ranges_list, areas=areas, stratified=stratified, maxtries=maxtries)

		# Stochastic mapping of events on the right branch
		events_table_for_branch_below_Right_node = stochastic_map_branch(nodenum_at_top_of_branch=right_branch_decnode, trtable=trtable, Qmat=Qmat, state_indices_0based=state_indices_0based, ranges_list=ranges_list, areas=areas, stratified=stratified, maxtries=maxtries)
	
		# Store the text representation
		# (extract to table with events_txt_into_events_table() )
		left_branch_events_txt = events_table_into_txt(events_table_for_branch_below_Left_node)
	
		right_branch_events_txt = events_table_into_txt(events_table_for_branch_below_Right_node)
	
		trtable$anagenetic_events_txt_below_node[left_branch_decnode] = left_branch_events_txt
		trtable$anagenetic_events_txt_below_node[right_branch_decnode] = right_branch_events_txt
		
# 		print(left_branch_decnode)
# 		print(left_branch_events_txt)
# 		print(right_branch_decnode)
# 		print(right_branch_events_txt)
	
		} # END for (j in edges_to_visit_uppass)
		  # (ENDING uppass loop)
	
	#print(trtable$anagenetic_events_txt_below_node)
	

	# Add cladogenetic events and re-arrange columns
	if (stratified == FALSE)
		{
		trtable = add_cladogenetic_events_to_trtable(trtable, BioGeoBEARS_run_object=res$inputs, tipranges, stratified=stratified, piecenum=NULL)
		} else {
		trtable = add_cladogenetic_events_to_trtable(trtable, BioGeoBEARS_run_object=res$inputs, tipranges, stratified=stratified, piecenum=piecenum)
		}

	#print(trtable$anagenetic_events_txt_below_node)

		
	# If stratified, don't rearrange
	if (stratified == FALSE)
		{
		first_colnums = 1:(ncol(trtable)-4)
		last4_colnums = (ncol(trtable)-3):(ncol(trtable))
		new_colnums = c(first_colnums, last4_colnums[c(2,3,4,1)])
		trtable = trtable[,new_colnums]
		}

	# Convert master_table_timeperiod_i to trtable
	if (stratified == TRUE)
		{
		#cat("\n\nprint(dim(trtable)):\n\n")
		#print(dim(trtable))

		#cat("\n\nprint(master_table_timeperiod_i[subtable_rownums,]):\n\n")
		#print(master_table_timeperiod_i[subtable_rownums,])

		
		# We may be modifying just PART of the subtable
		stochastic_mapping_inputs$master_table_timeperiod_i[subtable_rownums,] = trtable
		# Look at results
		stochastic_mapping_inputs$master_table_timeperiod_i
		
		return(stochastic_mapping_inputs$master_table_timeperiod_i)
		} else {
		# Look at results
		trtable[nodenums,]
		return(trtable)
		}
	
	return(stop("ERROR: you should not reach this."))
	} # END stochastic_map_given_inputs <- function(stochastic_mapping_inputs, piecenum=NULL, maxtries=40000, seedval=12345, include_null_range)

