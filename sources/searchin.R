##
# SEARCHIN Application Example
# ----------------------------------------------
# @author: Alessandro Vasciaveo
# @email: alessandro.vasciaveo@cumc.columbia.edu
# @copyright: 2019
# ----------------------------------------------

	isWorkspaceToClean <- TRUE
	if (isWorkspaceToClean)
		rm(list = ls())

	.isNeo4jServerOffline <- TRUE	
	set.seed(666)
	source("libs/searchin-class.R")	# Load SEARCHING libraries and utilities

	## Reading TXT file with the list of the ligands ----
  	tmp <- read_csv("input/putative_ligands_mouse_new_data_jul2019.txt",col_types = c("cccccccc"))
		ligands_list <- tmp$GeneSymbol
	## Reading CSV file with a table with the modualtors analysis ----
  	modulators_table <- read_csv("input/cindy-analysis.csv", col_types = "ccd")
	## Reading CSV file with a table with the VIPER analysis ----
  	receptors_table <- read_delim("input/viper-analysis.csv", delim = "\t" , col_types = "ccdddd")
  
	## Creating a SEARCHING object AND Running the analysis as pipeline ----
  my_sin <- SEARCHIN() %>% 
  	addLigandslist(ligands_list) %>% 
  	addModulatorTable(modulators_table) %>% 
  	addReceptorTable(receptors_table) %>% 
  	runModel() %>%
  	generateRank()
  
  my_sin
  
  require(pryr)
  object_size(my_sin)
  
  
  