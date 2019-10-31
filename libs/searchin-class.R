##
# SEARCHIN class definition
# -------------------------

rm(list = ls())

  require(methods)
  require(sloop)
  require(tidyverse)
  

  source("/Users/av2729/Workspace/SEARCHIN/libs/searchin-utility-class.R")
  
  source("/Users/av2729/Workspace/SEARCHIN/libs/preppiTable-class.R")
  source("/Users/av2729/Workspace/SEARCHIN/libs/viperTable-class.R")
  source("/Users/av2729/Workspace/SEARCHIN/libs/cindyTable-class.R")

  # setOldClass(c('CindyTable', 'PreppiTable', 'ViperTable'))  

  SEARCHIN = setClass( Class = "SEARCHIN",
    slots = c(
      preppi_obj = "PreppiTable" ,
      viper_obj = "ViperTable" ,
      cindy_obj = "CindyTable" ,
      
      ligand_list = "vector" ,
      
      id = "character" ,
      filename = "character" ,
      searchin_table = "tbl_df" ,
      private_table = "tbl_df"
    )
  )
  
  setMethod(
    f = "initialize", 
    signature = "SEARCHIN",
    function(.Object , ...) {
      .Object <- callNextMethod()
      # .Object@table <- tibble( rank = integer() , ligand = character() , receptor = character() , pvalue = double() )
      # .Object@table <- tibble( rank = integer() , interaction=character() , pvalue = double() )
      .Object@searchin_table <- tibble()
      .Object@private_table <- tibble()
      .Object
    }
  )
  
  setMethod(
    f = "show",
    signature = "SEARCHIN",
    definition = function(object) {
      cat(paste("An object of class", class(object), "\n"))
      cat("------------------------------\n")
      cat(">>> Using" , length(object@ligand_list) , "ligands\n")
      cat(">>> Using" , nrow(object@cindy_obj@cindy_table) , "modulators\n")
      cat(">>> Using" , nrow(object@viper_obj@viper_table) , "receptors\n")
      cat(">>> Using" , nrow(object@preppi_obj@preppi_table) , "PPIs interactions\n")
      # message(">>> ID: " , object@id)
      print(object@searchin_table)
      cat("------------------------------\n")
    }
  )
  
  setGeneric( "showAll", function(object) standardGeneric("showAll") )
  setMethod(
    f = "showAll",
    signature = "SEARCHIN",
    definition = function(object) {
      cat(">>> Showing all the tables componing SEARCHIN\n")
      cat("------------------------------\n")
      message(">>> Showing PrePPI Table:")
      print(object@preppi_obj)
      message(">>> Showing CINDy Table:")
      print(object@cindy_obj)
      message(">>> Showing VIPER Table:")
      print(object@viper_obj)      
      cat("------------------------------\n")
    }
  )  
  
  setValidity( Class = "SEARCHIN" ,
    validSearchinTable <- function(object) {
      if ( is_integer(object@table$rank) & is_character(object@table$ligand) & is_character(object@table$receptor) & is_double(object@table$pvalue) )
        TRUE
      else
        message("*** Error in SEARCHIN data format ***")
    })  
  
  setGeneric( "getTable", function(a) standardGeneric("getTable") )
  # setGeneric("setTable<-", function(object, value) standardGeneric("setTable<-") )
  
  setMethod("getTable", "SEARCHIN", function(a) a@table)

  # ## Cross-Species ID conversion
  # setGeneric( "crossSpeciesAnalysis", function(object,...) standardGeneric("crossSpeciesAnalysis") )
  # setMethod(
  #   f = "crossSpeciesAnalysis",
  #   signature = "SEARCHIN",
  #   definition = function(object,ligand_organism="rat",receptor_organism="mouse",report_table_organism="mouse") {
  #     cat(">>> Performing Cross-Species Analysis ...\n")
  #     if( is.null(object@preppi_obj) | is.null(object@viper_obj) | is.null(object@cindy_obj) )
  #     {
  #       message(">>> Missing Tables. Load tables first.")
  #       return(0)
  #     }  
  #     
  #     object@preppi_obj
  #     object@uniprot_mapping
  #     
  #     return(object)
  #   }
  # )
  
  ## Add Ligand List
  setGeneric( "addLigandslist", function(object,...) standardGeneric("addLigandslist") )
  setMethod(
    f = "addLigandslist",
    signature = "SEARCHIN",
    definition = function(object,ligand_list,ligand_organism="mouse",ligand_id_type="entrezid") {
      
      cat(">>> >>> Processed" , length(ligand_list) , "ligands\n")
      object@ligand_list <- ligand_list
      cat(">>> Adding List of Ligands ...\n") 
      # slot(object,"preppi_obj") <<- object@preppi_obj
      
      invisible(object)
    }
  )   
  
  ## Add Modulator Table
  setGeneric( "addModulatorTable", function(object,...) standardGeneric("addModulatorTable") )
  setMethod(
    f = "addModulatorTable",
    signature = "SEARCHIN",
    definition = function(object,modulators_table,receptor_organism="mouse",receptor_id_type="entrezid") {
      
      cat(">>> >>> Processed" , nrow(modulators_table) , "modulators\n")
      setCindyTable(object@cindy_obj) <- modulators_table
      cat(">>> Adding Modulator Table ...\n") 
      # slot(object,"preppi_obj") <<- object@preppi_obj
      
      invisible(object)
    }
  ) 
  
  ## Add Receptor Table
  setGeneric( "addReceptorTable", function(object,...) standardGeneric("addReceptorTable") )
  setMethod(
    f = "addReceptorTable",
    signature = "SEARCHIN",
    definition = function(object,receptor_table,receptor_organism="mouse",receptor_id_type="entrezid") {
      
      cat(">>> >>> Processed" , nrow(receptor_table) , "receptors\n")
      setViperTable(object@viper_obj) <- receptor_table
      cat(">>> Adding Receptor Table ...\n") 
      # slot(object,"preppi_obj") <<- object@preppi_obj
      
      invisible(object)
    }
  )    
  
  ## Running SEARCHING Algorithm
  setGeneric( "runModel", function(object,...) standardGeneric("runModel") )
  setMethod(
    f = "runModel",
    signature = "SEARCHIN",
    definition = function(object,ligand_organism="rat",receptor_organism="mouse",report_table_organism="mouse") {
      cat(">>> Performing Cross-Species Analysis ...\n")
      if( is.null(object@preppi_obj) | is.null(object@viper_obj) | is.null(object@cindy_obj) )
      {
        message(">>> Missing Tables. Load tables first.")
        return(0)
      }  
      
      s <- SUtils()
      
      l <- convertFeatures(object@ligand_list,"SYMBOL","ENSEMBL","mouse")
      l <- s %>% convert_ens( what = l, from = "mouse" , to = "human" ) # Here is the Cross Species conversion. Change the name of the function
      l <- s@uniprot_mapping$uniprotswissprot[ match( l , s@uniprot_mapping$ensembl_gene_id) ] # Get Protein ID (Uniprot) from Human ENSEMBL ID
      
      tmp <- getViperTable(object@viper_obj)
      r <- convertFeatures( tmp$Symbol ,"SYMBOL","ENSEMBL","mouse")
      r <- s %>% convert_ens( what = r, from = "mouse" , to = "human" ) # Here is the Cross Species conversion. Change the name of the function
      r <- s@uniprot_mapping$uniprotswissprot[ match( r , s@uniprot_mapping$ensembl_gene_id) ] # Get Protein ID (Uniprot) from Human ENSEMBL ID
      
      # s@uniprot_mapping[ s@uniprot_mapping$uniprotswissprot %in% "P14174", ]
      
      tmp <- runPreppiModel( l , r )
      setPreppiTable(object@preppi_obj) <- tmp
      
      object@private_table <- getPreppiTable(object@preppi_obj)
      object@private_table$p_ligand_sym <- s@uniprot_mapping[ match( object@private_table$p1 , s@uniprot_mapping$uniprotswissprot ) , "symbol" , drop = TRUE ]
      object@private_table$p_receptor_sym <- s@uniprot_mapping[ match( object@private_table$p2 , s@uniprot_mapping$uniprotswissprot ) , "symbol" , drop = TRUE ]
      
      x <- object@private_table$p_ligand_sym
      x <- x %>% convertFeatures("SYMBOL","ENSEMBL","human")
      x <- s %>% convert_ens( what = x , from = "human" , to = "mouse" )
      x <- x %>% convertFeatures("ENSEMBL","SYMBOL","mouse")
      object@private_table$p_ligand_sym_mouse <- x
      
      x <- object@private_table$p_receptor_sym
      x <- x %>% convertFeatures("SYMBOL","ENSEMBL","human")
      x <- s %>% convert_ens( what = x , from = "human" , to = "mouse" )
      x <- x %>% convertFeatures("ENSEMBL","SYMBOL","mouse")
      object@private_table$p_receptor_sym_mouse <- x
      
      invisible(object)
    }
  )      
  
  ## Cross-Species ID conversion
  setGeneric( "generateRank", function(object,...) standardGeneric("generateRank") )
  setMethod(
    f = "generateRank",
    signature = "SEARCHIN",
    definition = function(object,ligand_organism="rat",receptor_organism="mouse",report_table_organism="mouse") {
      require(RobustRankAggreg)
      cat("- Integrating Evidences - Using Rank Aggregation (RRA) \n" )
      if( is.null(object@private_table) )
      {
        message(">>> Missing PPIs Tables. Run PreppiModel first.\n")
        return(0)
      }  
      
      .table <- object@private_table
      
      t_1 <- .table %>% dplyr::rename(preppi_pvalue=pvalue)
      t_2 <- getCindyTable(object@cindy_obj) %>% dplyr::rename(cindy_pvalue=pvalue)
      .table <- inner_join( t_1 , t_2 , by = c("p_receptor_sym_mouse"="Modulator") )
      
      t_3 <- getViperTable(object@viper_obj) %>% dplyr::rename(viper_pvalue=pvalue)
      .table <- inner_join( .table , t_3, by = c("p_receptor_sym_mouse"="Symbol") )
      
      object@private_table <- .table
      
      cat(">>> Using order statistics as meta-analysis technique\n")
      
      preppi_order <- rev( order( .table$preppi_pvalue , decreasing = TRUE ) )
      a <- paste( .table$p_ligand_sym_mouse[ preppi_order ] , .table$p_receptor_sym_mouse[ preppi_order ] , sep = "-" )
      
      viper_order <- rev( order( .table$viper_pvalue , decreasing = TRUE ) )
      b <- paste( .table$p_ligand_sym_mouse[ viper_order ] , .table$p_receptor_sym_mouse[ viper_order ] , sep = "-" )
      
      cindy_order <- rev( order( .table$cindy_pvalue , decreasing = TRUE ) )
      c <- paste( .table$p_ligand_sym_mouse[ cindy_order ] , .table$p_receptor_sym_mouse[ cindy_order ] , sep = "-" )
      
      glist <- list( a , b , c )
      
      ranked_interactions <- aggregateRanks( glist = glist, N = mean( sapply( glist , length ) ) , method = "RRA" )
      rownames(ranked_interactions) <- NULL
      object@searchin_table <- as_tibble(ranked_interactions) %>% 
        dplyr::rename(interaction=Name,pvalue=Score) %>% 
        mutate(rank=1:n())
      # tibble( rank = integer() , ligand = character() , receptor = character() , pvalue = double() )
      # print(ranked_interactions[1:10,])
      
      invisible(object)
    }
  )

##
# Loading Data
# ------------

#   message(">>> Reading data from a Excel format ...")
#   # detach("package:xlsxjars",unload = TRUE)
#   # detach("package:xlsx",unload = TRUE)
#   # detach("package:rJava",unload = TRUE)
#   options( java.parameters = "-Xms4096m" )
#   options( java.parameters = "-Xmx8000m" )
#   library(rJava)
#   library(xlsx)
#   library(tidyverse)
# 
# 		cindy_filename <- "data/input/excel/cindy-pvalues.xlsx"
#     cindy_dt <- as_tibble( read.xlsx2( file = cindy_filename , sheetIndex = 1 ) )
#     cindy_dt$p.value <- as.numeric( as.character( cindy_dt$p.value ) )
#     cindy_dt$Modulator <- as.character( cindy_dt$Modulator )
#     setnames(cindy_dt,old = "p.value", new = "cindy.p.value" )
#     setkey(cindy_dt,Modulator)
#     
# 		viper_filename <- "data/input/excel/viper-pvalues.xlsx"
#     viper_dt <- as.data.table(read.xlsx2( file = viper_filename , sheetIndex = 1 ))
#     viper_dt$p.value <- as.numeric( as.character( viper_dt$p.value ) )
#     setnames(viper_dt,old = "p.value", new = "viper.p.value" )
#     setkey(viper_dt,Receptor)
#     
# 		preppi_filename <- "data/input/excel/preppi-pvalues.xlsx"
#     preppi_dt <- as.data.table(read.xlsx2( file = preppi_filename , sheetIndex = 1 ))
#     preppi_dt$p.value <- as.numeric( as.character( preppi_dt$p.value ) )
#     setnames(preppi_dt,old = "p.value", new = "preppi.p.value" )
#     setkey(preppi_dt,Ligand.Gene.Name)
  
    #   colnames(ranked_interactions)[1] <- "Ligand-Receptor Predicted Interaction"
    # 
    # filename <- "data/output/excel/searchin-table.xlsx"
    # message("- Saving SEARCHIN Table in: " , filename )
    # write.xlsx( x = ranked_interactions , file = filename , sheetName = "SEARCHIN Algorithm Output" , row.names = FALSE )
  
  
  
  
  