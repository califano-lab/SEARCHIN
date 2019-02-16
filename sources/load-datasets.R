##
# Loading Data
# ------------

  message(">>> Reading data from a Excel format ...")
  # detach("package:xlsxjars",unload = TRUE)
  # detach("package:xlsx",unload = TRUE)
  # detach("package:rJava",unload = TRUE)
  options( java.parameters = "-Xms4096m" )
  options( java.parameters = "-Xmx8000m" )
  library(rJava)
  library(xlsx)
  library(data.table)

		cindy_filename <- "data/input/excel/cindy-pvalues.xlsx"
    cindy_dt <- as.data.table( read.xlsx2( file = cindy_filename , sheetIndex = 1 ) )
    cindy_dt$p.value <- as.numeric( as.character( cindy_dt$p.value ) )
    setnames(cindy_dt,old = "p.value", new = "cindy.p.value" )
    setkey(cindy_dt,Modulator)
    
		viper_filename <- "data/input/excel/viper-pvalues.xlsx"
    viper_dt <- as.data.table(read.xlsx2( file = viper_filename , sheetIndex = 1 ))
    viper_dt$p.value <- as.numeric( as.character( viper_dt$p.value ) )
    setnames(viper_dt,old = "p.value", new = "viper.p.value" )
    setkey(viper_dt,Receptor)
    
		preppi_filename <- "data/input/excel/preppi-pvalues.xlsx"
    preppi_dt <- as.data.table(read.xlsx2( file = preppi_filename , sheetIndex = 1 ))
    preppi_dt$p.value <- as.numeric( as.character( preppi_dt$p.value ) )
    setnames(preppi_dt,old = "p.value", new = "preppi.p.value" )
    setkey(preppi_dt,Ligand.Gene.Name)
    
		