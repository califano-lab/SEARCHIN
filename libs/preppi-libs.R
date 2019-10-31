##
# Functions for Ligand/Receptor interaction estimation
# ----------------------------------------------------
# source("/Users/av2729/Workspace/SEARCHIN/libs/preppi-libs.R")

# require(RNeo4j)
require(qvalue)
require(data.table)
# require(tidyverse)

setClass("PrePPIClass",
  slots = list( name = "character" , age = "numeric" , preppi_dt = "data.table" , cache_table = "tbl" ) )

  ##
  # Functions to store/load PPI pvalues in caches
  # ---------------------------------------------
  load_ppi_cache <- function(path)
  {
    x <- read_csv(path,col_types = c("ccd"))
    # ppi_cache <- tibble( p1 = character() , p2 = character() , pvalue = double() )
    return( x )
  }
  store_ppi_cache <- function(x,path)
  {
    return( write_csv(x,path) )
  }
    
  get_ppi_pvalue <- function(cache_table,arg0,arg1) { 
    res <- cache_table %>% filter( p1 == arg0 & p2 == arg1) %>% distinct() %>% dplyr::select("pvalue") %>% pull()
    res <- as.double(res)
    return(res)
  }
  
  ##
  # Function to generate a protien family specific pvalue for PrePPI
  # ----------------------------------------------------------------
  # debug: works only when cache_table is NULL
  # cache_table: is not a paramter but a object with global scope
  message("- Loading Function to Generate PrePPI Null Model through GraphDB PPI Interactome ")
    getPrePPIScoresFromGDB <- function( p1 , p2 , preppi_obj , debug = FALSE )
    {
      pvalue_to_return <- get_ppi_pvalue( preppi_obj@cache_table , p1 , p2 )
      
      if ( !is.na(pvalue_to_return) & is.double(pvalue_to_return) ) # the value is a number (double): Then is consistent to the format of a pvalue
      {
        return(pvalue_to_return)
        
      } else { # ppi is not present in cache we need to perform the query to GraphDB, build the null model and generate the pvalue
        
        # query.p1 <- paste0( "MATCH (p1:PreppiProtein)-[r:INTERACTS_WITH]->(p2:PreppiProtein) WHERE p1.name={p1} RETURN p1.name,r.final_score,p2.name" )
        # result.p1 <- cypher(graph, query.p1 , p1=p1)
        # query.p2 <- paste0( "MATCH (p1:PreppiProtein)-[r:INTERACTS_WITH]->(p2:PreppiProtein) WHERE p2.name={p2} RETURN p1.name,r.final_score,p2.name" )
        # result.p2 <- cypher(graph, query.p2 , p2=p2)
        query.p1 <- paste0( "MATCH (p1:Protein)-[r:interacts_with]->(p2:Protein) WHERE p1.name={p1} RETURN p1.name,r.final_score,p2.name" )
        result.p1 <- cypher(graph, query.p1 , p1=p1)
        query.p2 <- paste0( "MATCH (p1:Protein)-[r:interacts_with]->(p2:Protein) WHERE p2.name={p2} RETURN p1.name,r.final_score,p2.name" )
        result.p2 <- cypher(graph, query.p2 , p2=p2)
        
        stopifnot( identical( result.p1[ result.p1$p2.name %in% p2 , "r.final_score" ] , result.p2[ result.p2$p1.name %in% p1 , "r.final_score" ] ) )
        
        .test <- result.p1[ result.p1$p2.name %in% p2 , "r.final_score" ]
        
        x <- result.p1$r.final_score
        y <- result.p2$r.final_score
        z <- rep(0,2*20259-length(x)-length(y))
        print( paste0("- Null model for proteins " , p1 , " -> " , p2 , " using interactions p1:" , length(x) , " , p2:" , length(y) , " , miss:" , length(z) , " , tot:" , sum(length(x),length(y),length(z))) )
        
        x <- x/(x+600)
        y <- y/(y+600)
        .test <- .test/(.test+600)
        
        pvalue_to_return <- empPvals( .test , c(x,y,z) )
        
        if ( debug )
          return( list( p.value = pvalue_to_return , p.test = .test , p1.values = x , p2.values = y ) )
        else
          return( pvalue_to_return )
      } 
    }  