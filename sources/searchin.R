##
# SEARCHIN Application Example
# ----------------------------

	isWorkspaceToClean <- TRUE
	if (isWorkspaceToClean)
		rm(list = ls())
	
	set.seed(666)
	source("libs/searchin-class.R")	

	## Reading TXT file with the list of the ligands ----
  tmp <- read_csv("/Users/av2729/Workspace/SEARCHIN/input/putative_ligands_mouse_new_data_jul2019.txt",col_types = c("cccccccc"))
	ligands_list <- tmp$GeneSymbol
	## Reading CSV file with a table with the modualtors analysis ----
  modulators_table <- read_csv("/Users/av2729/Workspace/SEARCHIN/input/cindy-analysis.csv", col_types = "ccd")
	## Reading CSV file with a table with the VIPER analysis ----
  receptors_table <- read_delim("/Users/av2729/Workspace/SEARCHIN/input/viper-analysis.csv", delim = "\t" , col_types = "ccdddd")
  
	## Creating a SEARCHING object ----
  my_sin <- SEARCHIN()
	
  my_sin <- my_sin %>% 
  	addLigandslist(ligands_list) %>% 
  	addModulatorTable(modulators_table) %>% 
  	addReceptorTable(receptors_table) %>% 
  	runModel() %>%
  	generateRank()
  
  my_sin
  
  require(pryr)
  object_size(my_sin)
  
  
  