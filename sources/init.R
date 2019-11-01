##
# Environment Initialization
# ----------------------------------------------
# @author: Alessandro Vasciaveo
# @email: alessandro.vasciaveo@cumc.columbia.edu
# @copyright: 2019
# ----------------------------------------------


# From CRAN Repository
  list.of.required.packages <- c( "data.table" , "RobustRankAggreg" )
  list.of.new.packages <- list.of.required.packages[ !(list.of.required.packages %in% installed.packages()[,"Package"]) ]
  if( length(list.of.new.packages) > 0 ) install.packages(list.of.new.packages)
  