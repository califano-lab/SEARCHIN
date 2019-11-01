##
# ViperTable class definition
# ----------------------------------------------
# @author: Alessandro Vasciaveo
# @email: alessandro.vasciaveo@cumc.columbia.edu
# @copyright: 2019
# ----------------------------------------------

  require(methods)
  require(sloop)

  setOldClass(c('tbl_df', 'tbl', 'data.frame','tibble'))
  
  ViperTable = setClass( Class = "ViperTable",
    slots = c(
      id = "character" ,
      filename = "character" ,
      viper_table = "tbl_df"
    )
  )
  
  setMethod(
    f = "initialize", 
    signature = "ViperTable",
    function(.Object , ...) {
      .Object <- callNextMethod()
      .Object@viper_table <- tibble(GeneID = character() , Symbol = character() , NES = double() , pvalue = double() ,
        FDR = double() , Bonferroni = double() )
      # if (is.null(id))
      #   .Object@id = timestamp(quiet = T)
      # else
      #   .Object@id = id      
      .Object
    }
  )
  
  setMethod(
    f = "show",
    signature = "ViperTable",
    definition = function(object) {
      cat(paste("An object of class", class(object), "\n"))
      cat("------------------------------\n")
      message(">>> ID: " , object@id)
      print(object@viper_table)
      cat("------------------------------\n")
    }
  )
  
  setValidity( Class = "ViperTable" ,
    validViperTable <- function(object) {
      if ( is_character(object@viper_table$GeneID) & is_character(object@viper_table$Symbol) & is_double(object@viper_table$NES) & is_double(object@viper_table$pvalue) )
        TRUE
      else
        message("*** Error in ViperTable data format ***")
    })  
  
  setGeneric( "getViperTable", function(object) standardGeneric("getViperTable") )
  setGeneric("setViperTable<-", function(object, value) standardGeneric("setViperTable<-") )
  
  setMethod("getViperTable", "ViperTable", function(object) object@viper_table)
  setMethod("setViperTable<-", "ViperTable", function(object, value) {
    object@viper_table <- value
    validViperTable(object)
    object
  })  
