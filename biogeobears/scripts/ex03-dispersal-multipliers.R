#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/biogeogbears-utilities.R")
load.biogeobears()

###################################################################
## Set up the phylogeny and geography
phylogeny.path = "data/Psychotria_5.2.newick"
psychotria.tree = read.tree(phylogeny.path)
geography.path = "data/Psychotria_geog.data"
dec.ss.run = define_BioGeoBEARS_run()
dec.ss.run$trfn = phylogeny.path
dec.ss.run$geogfn = geography.path


###################################################################
## Set up the dispersal restrictions
dec.ss.run$dispersal_multipliers_fn = "data/Psychotria-steppingstone-dispersal.txt"
# This function loads the dispersal multiplier matrix etc. from the text files
# into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
dec.ss.run = readfiles_BioGeoBEARS_run(dec.ss.run)

# Configure the run
dec.ss.run = configure.standard.biogeobears.run(dec.ss.run)

## Check
check_BioGeoBEARS_run(dec.ss.run)

## Execute the analysis
dec.results = bears_optim_run(dec.ss.run)
save(dec.results, file="results/Psychotria.DEC.steppingstone.Rdata")
