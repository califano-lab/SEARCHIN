##
# Data Integration
# ----------------

library(RobustRankAggreg)

	x <- merge( cindy_dt , preppi_dt , by.x = "Modulator" , by.y = "Receptor.Gene.Name" )
	searchin.table <- merge( x , viper_dt , by.x = "Modulator" , by.y = "Receptor" )

  message("- Using order statistics as meta-analysis technique")
  {
    
    preppi.order <- rev( order( searchin.table$preppi.p.value , decreasing = TRUE ) )
    a <- paste( searchin.table$Ligand.Gene.Name[ preppi.order ] , searchin.table$Modulator[ preppi.order ] , sep = "-" )

    viper.order <- rev( order( searchin.table$viper.p.value , decreasing = TRUE ) )
    b <- paste( searchin.table$Ligand.Gene.Name[ viper.order ] , searchin.table$Modulator[ viper.order ] , sep = "-" )

    cindy.order <- rev( order( searchin.table$cindy.p.value , decreasing = TRUE ) )
    c <- paste( searchin.table$Ligand.Gene.Name[ cindy.order ] , searchin.table$Modulator[ cindy.order ] , sep = "-" )
    
    glist <- list( a , b , c )
    
    message("- Integrating Evidences - Rank Aggregation (RRA) " )
    ranked_interactions <- aggregateRanks( glist = glist, N = mean( sapply( glist , length ) ) , method = "RRA" )
    # View(ranked_interactions)
    
    colnames(ranked_interactions)[1] <- "Ligand-Receptor Predicted Interaction"
    
    filename <- "data/output/excel/searchin-table.xlsx"
    message("- Saving SEARCHIN Table in: " , filename )
    write.xlsx( x = ranked_interactions , file = filename , sheetName = "SEARCHIN Algorithm Output" , row.names = FALSE )

  }
  