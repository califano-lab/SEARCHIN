# elapsed.time <- system.time( source("sources/analysis-for-paper/searchin.R") )
  
	isWorkspaceToClean <- TRUE
	if (isWorkspaceToClean)
		rm(list = ls())
  
  set.seed(666)
  # source("../vaxtools/R/utils.R")
  # source("../vaxtools/R/cross-species-utils.R")
  # entrez2geneSymbol.mapped <- entrez2geneSymbol.load( organism = "mm" )

##
# SEARCHIN Algorithm
# ------------------

## Load Datasets
  source("sources/load-datasets.R")
  source("sources/evidence-integrator.R")
  
  message("\n----------------------------------------  *** Analysis COMPLETE *** ---------------------------------------- ")
    # message("Total time: " , elapsed.time[3]/60 , " minutes" )    