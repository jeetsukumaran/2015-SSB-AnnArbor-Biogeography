#!/usr/bin/env Rscript

## Load BioGeoBEARS and apply patches
source("scripts/biogeogbears-utilities.R")
load.biogeobears()

## Create the run object for the DECJ+J analysis
decj.run = define_BioGeoBEARS_run()

## Inspect the phylogeny
phylogeny.path = "data/Psychotria_5.2.newick"
moref(phylogeny.path)
psychotria.tree = read.tree(phylogeny.path)
plot(psychotria.tree)
prt(psychotria.tree, get_tipnames=TRUE)

## Set the phylogeny
decj.run$trfn = phylogeny.path
decj.run$trfn

## Inspect the geography
geography.path = "data/Psychotria_geog.data"
moref(geography.path)
getranges_from_LagrangePHYLIP(lgdata_fn=geography.path)

## Set the geography
decj.run$geogfn = geography.path
decj.run$geogfn

## Add the founder-event dispersal parameter
decj.run$BioGeoBEARS_model_object@params_table["j","type"] = "free"
decj.run$BioGeoBEARS_model_object@params_table["j","init"] = 0.0001
decj.run$BioGeoBEARS_model_object@params_table["j","est"] = 0.0001

## Set a standard run configuration
## Note, due to the way R passes objects, you need to
## rebind the name "decj.run" to the returned object
decj.run = configure.standard.biogeobears.run(decj.run)

## Check
check_BioGeoBEARS_run(decj.run)

## Execute the analysis
decj.results = bears_optim_run(decj.run)
save(decj.results, file="results/Psychotria.DECJ.Rdata")

## View the results by ranges
plot.biogeobears.results.ranges(
        results.object=decj.results,
        plot.type="text",
        analysis.title="DECJ"
        )

## View the results by areas
plot.biogeobears.results.areas(
        results.object=decj.results,
        analysis.title="DECJ"
        )

## Get the results as a table
decj.range.table = get.biogeobears.results.by.range.table(decj.results)
write.table(decj.range.table, "results/Psychotria.DECJ.ranges.tsv", sep="\t", row.names=F)

## Get the area results as a table
decj.area.table = get.biogeobears.results.by.area.table(decj.results)
write.table(decj.area.table, "results/Psychotria.DECJ.areas.tsv", sep="\t", row.names=F)

## Look up a node
ndi = getMRCA(psychotria.tree, c("P_hawaiiensis_Makaopuhi","P_wawraeDL7428"))
decj.range.table[ndi,]
