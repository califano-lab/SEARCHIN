##
# PreppiTable class definition
# ----------------------------------------------
# @author: Alessandro Vasciaveo
# @email: alessandro.vasciaveo@cumc.columbia.edu
# @copyright: 2019
# ----------------------------------------------

	source("libs/ppi-null-model-generator.R")

  require(methods)
  require(sloop)

  setOldClass(c('tbl_df', 'tbl', 'data.frame','tibble'))
  
  PreppiTable = setClass( Class = "PreppiTable",
    slots = c(
      id = "character" ,
      filename = "character" ,
      preppi_table = "tbl_df"
    )
  )
  
  setMethod(
    f = "initialize", 
    signature = "PreppiTable",
    function(.Object , ...) {
      .Object <- callNextMethod()
      .Object@preppi_table <- tibble( p1 = character() , p2 = character() , pvalue = double() )
      # if (is.null(id))
      #   .Object@id = timestamp(quiet = T)
      # else
      #   .Object@id = id      
      .Object
    }
  )
  
  setMethod(
    f = "show",
    signature = "PreppiTable",
    definition = function(object) {
      cat(paste("An object of class", class(object), "\n"))
      cat("------------------------------\n")
      message(">>> ID: " , object@id)
      print(object@preppi_table)
      cat("------------------------------\n")
    }
  )
  
  setValidity( Class = "PreppiTable" ,
    validPreppiTable <- function(object) {
      if ( is_character(object@preppi_table$p1) & is_character(object@preppi_table$p2) & is_double(object@preppi_table$pvalue) )
        TRUE
      else
        message("*** Error in PreppiTable data format ***")
    })  
  
  setGeneric( "getPreppiTable", function(object) standardGeneric("getPreppiTable") )
  setGeneric("setPreppiTable<-", function(object, value) standardGeneric("setPreppiTable<-") )
  
  setMethod("getPreppiTable", "PreppiTable", function(object) object@preppi_table)
  setMethod("setPreppiTable<-", "PreppiTable", function(object, value ) {
    object@preppi_table <- value
    validPreppiTable(object)
    object
  })  
  