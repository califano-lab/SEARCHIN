##
# SEARCHIN Utility class definition
# ---------------------------------
# source("/Users/av2729/Workspace/SEARCHIN/libs/searchin-utility-class.R")

  SUtils = setClass( Class = "SUtils",
    slots = c(
      id = "character" ,
      filename = "character" ,
      table = "tbl_df" ,
      uniprot_mapping = "tbl_df" ,
      ortho_mouse2human = "tbl_df" ,
      ortho_human2mouse = "tbl_df" 
      # entrez2sym = "tbl_df"
      # ortho_rat2mouse = "tbl_df" 
    )
  )

  setMethod(
    f = "initialize", 
    signature = "SUtils",
    function(.Object , ...) {
      .Object <- callNextMethod()
      # .Object@table <- tibble( rank = integer() , ligand = character() , receptor = character() , pvalue = double() )
      .Object@uniprot_mapping <- read_csv("data/uniprot-entrezgeneid-symbol-mapping-human-mouse-rat.csv" , col_types = c("ccccccc"))
      .Object@ortho_mouse2human <- read_csv( "data/mouse2human.csv" , col_types = c("ccc"))
      .Object@ortho_human2mouse <- read_csv( "data/human2mouse.csv" , col_types = c("ccc"))
      # .Object@entrez2sym <- read_csv("/Users/av2729/Workspace/vaxtools/data/mapping-human-mouse-entrez-id-gene-symbols-gene-biotype.csv" , col_types = c("cccc"))
      # .Object@ortho_rat2mouse <- read_csv( "/Users/av2729/Workspace/SEARCHIN/utils/rat2mouse.csv" , col_types = c("ccc"))
      # .Object@ortho_rat2mouse_entrez <- read_csv( "/Users/av2729/Workspace/SEARCHIN/utils/rat2mouse-entrez.csv" , col_types = c("ccc"))
      .Object
    }
  )  
  
  setMethod(
    f = "show",
    signature = "SUtils",
    definition = function(object) {
      cat(paste("An object of class", class(object), "\n"))
      cat("------------------------------\n")
      cat("uniprot table contains" , nrow(object@uniprot_mapping) , "\n")
      cat("ortho_mouse2human table contains" , nrow(object@ortho_mouse2human) , "\n")
      cat("ortho_human2mouse table contains" , nrow(object@ortho_human2mouse) , "\n")
      # cat("ortho_rat2mouse table contains" , nrow(object@ortho_rat2mouse) , "\n")
      # message(">>> ID: " , object@id)
      # print(object@uniprot_mapping)
      cat("------------------------------\n")
    }
  ) 
  
  # setGeneric( "convert_sym", function(object,what,from,to) standardGeneric("convert_sym") )
  # setMethod(
  #   f = "convert_sym",
  #   signature = c("SUtils") ,
  #   definition = function(object,what,from,to) {
  #     # what_converted <- object@uniprot_mapping$symbol[ match( object@preppi_obj@preppi_table$p1 ,
  #     #   object@uniprot_mapping$uniprotswissprot ) ]
  #     from <- base::match.arg(from, c("human","mouse"))
  #     to <- base::match.arg(to, c("human","mouse"))
  #     db_list <- c("mouse2human","human2mouse")
  #     db_name <- paste0(from,"2",to)
  #     if ( !(db_name %in% db_list) )
  #     {
  #       print(">>> convert_ens : *** Error *** : Conversion table not present")
  #       return(object)
  #     }
  # 
  #     db_selected <- switch( db_name ,
  #       "mouse2human" = object@ortho_mouse2human ,
  #       "human2mouse" = object@ortho_human2mouse )
  # 
  #     what_converted <- db_selected[ match( what , db_selected[,1] %>% pull() ) , 2 ] %>% pull()
  # 
  #     what_converted
  #   }
  # )
  
  setGeneric( "convert_ens", function(object,what,from,to) standardGeneric("convert_ens") )
  setMethod(
    f = "convert_ens",
    signature = c("SUtils") ,
    definition = function(object,what,from,to) {
      # what_converted <- object@uniprot_mapping$symbol[ match( object@preppi_obj@preppi_table$p1 ,
      #   object@uniprot_mapping$uniprotswissprot ) ]
      from <- base::match.arg(from, c("human","mouse"))
      to <- base::match.arg(to, c("human","mouse"))
      db_list <- c("mouse2human","human2mouse")
      db_name <- paste0(from,"2",to)
      if ( !(db_name %in% db_list) )
      {
        print(">>> convert_ens : *** Error *** : Conversion table not present")
        return(object)
      }

      db_selected <- switch( db_name ,
        "mouse2human" = object@ortho_mouse2human ,
        "human2mouse" = object@ortho_human2mouse )

      what_converted <- db_selected[ match( what , db_selected[,1] %>% pull() ) , 2 ] %>% pull()

      what_converted
    }
  )
  
  # s <- SUtils()
  # s
  # x <- read_csv( "/Users/av2729/Workspace/SEARCHIN/utils/mouse2human.csv" , col_types = c("ccc"))[1:10,1]
  # s %>% convert_ens( what = x$mouse_ens_id , from = "mouse" , to = "human" )
  
  require(org.Hs.eg.db)
  require(org.Mm.eg.db)
  convertFeatures <- function( ids , from = "ENSEMBL" , to = "SYMBOL" , organism = "human" )
  {
  	db = if( organism == "human" ) org.Hs.eg.db else org.Mm.eg.db 
  	suppressMessages( suppressWarnings( mapping <- mapIds( db , keys = ids , column = to , keytype = from , multiVals = "first" ) ) )
  	x <- mapping[ match( ids , names(mapping) ) ]
  	x <- as.character(x)
  	x <- ifelse( grepl("NULL",x) | is.na(x) , NA , x )
  	return( x )
  }
  
  
  
  
  
  
  
  
  