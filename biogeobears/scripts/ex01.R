#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/biogeogbears-utilities.R")
load.biogeobears()

## Create the run object for the DEC analysis
dec.run = define_BioGeoBEARS_run()

## Inspect the phylogeny
phylogeny.path = "data/Psychotria_5.2.newick"
moref(phylogeny.path)
psychotria.tree = read.tree(phylogeny.path)
plot(psychotria.tree)
prt(psychotria.tree, get_tipnames=TRUE)

## Set the phylogeny
dec.run$trfn = phylogeny.path
dec.run$trfn

## Inspect the geography
geography.path = "data/Psychotria_geog.data"
moref(geography.path)
getranges_from_LagrangePHYLIP(lgdata_fn=geography.path)

## Set the geography
dec.run$geogfn = geography.path
dec.run$geogfn

## Set a standard run configuration
## Note, due to the way R passes objects, you need to
## rebind the name "dec.run" to the returned object
dec.run = configure.standard.biogeobears.run(dec.run)

## Execute the analysis
dec.results = bears_optim_run(dec.run)
save(dec.results, file="results/Psychotria.DEC.Rdata")

## View the results
plot.biogeobears.results(
        results.object=dec.results,
        plot.type="text",
        analysis.title="DEC"
        )

## Get the results as a table
dec.table = get.biogeobears.results.table(dec.results)
write.table(dec.table, "results/Psychotria.DEC.tsv", sep="\t", row.names=F)

## Look up a node
ndi = getMRCA(psychotria.tree, c("P_hawaiiensis_Makaopuhi","P_wawraeDL7428"))
dec.table[ndi,]
