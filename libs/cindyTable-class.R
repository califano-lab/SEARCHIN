##
# CindyTable class definition
# ----------------------------------------------
# @author: Alessandro Vasciaveo
# @email: alessandro.vasciaveo@cumc.columbia.edu
# @copyright: 2019
# ----------------------------------------------

  require(methods)
  require(sloop)

  setOldClass(c('tbl_df', 'tbl', 'data.frame','tibble'))
  
  CindyTable = setClass( Class = "CindyTable",
    slots = c(
      id = "character" ,
      filename = "character" ,
      cindy_table = "tbl_df"
    )
  )
  
  setMethod(
    f = "initialize", 
    signature = "CindyTable",
    function(.Object , ...) {
      .Object <- callNextMethod()
      .Object@cindy_table <- tibble(Modulator = character() , TF = character() , 
      	# significantTriplets = integer() , 
      	pvalue = double())
      # if (is.null(id))
      #   .Object@id = timestamp(quiet = T)
      # else
      #   .Object@id = id      
      .Object
    }
  )
  
  setMethod(
    f = "show",
    signature = "CindyTable",
    definition = function(object) {
      cat(paste("An object of class", class(object), "\n"))
      cat("------------------------------\n")
      message(">>> ID: " , object@id)
      print(object@cindy_table)
      cat("------------------------------\n")
    }
  )
  
  setValidity( Class = "CindyTable" ,
    validCindyTable <- function(object) {
      if ( is_character(object@cindy_table$Modulator) & is_character(object@cindy_table$TF) & 
      		# is_integer(object@cindy_table$significantTriplets) & 
      		is_double(object@cindy_table$pvalue) )
        TRUE
      else
        message("*** Error in CindyTable data format ***")
    })  
  
  setGeneric( "getCindyTable", function(object) standardGeneric("getCindyTable") )
  setGeneric("setCindyTable<-", function(object, value) standardGeneric("setCindyTable<-") )
  
  setMethod("getCindyTable", "CindyTable", function(object) object@cindy_table)
  setMethod("setCindyTable<-", "CindyTable", function(object, value) {
    object@cindy_table <- value
    validCindyTable(object)
    object
  })  
  