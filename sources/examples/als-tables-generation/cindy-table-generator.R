##
# CINDy Null Model Generation
## ---------------------------------------------
# @author: Alessandro Vasciaveo
# @email: alessandro.vasciaveo@cumc.columbia.edu
# @copyright: 2019
# ----------------------------------------------

library(data.table)
library(dplyr)

  message("- Reading Nfkb1 Triplets to generate Null Model")
  cindy_analysis.env <- list()
  cindy_analysis.env$triplets.nfkb1 <- fread( "data/examples/als-tables-generation/triplets_Nfkb1.txt" , header = TRUE , stringsAsFactors = FALSE , strip.white = TRUE )

  cindy_analysis.env$triplets.nfkb1[ , n_triplets := length(Target) , by = list(Modulator) ]
  cindy_analysis.env$triplets.nfkb1[ , CMI_avg := mean(CMI) , by = list(Modulator) ]
  
  message(" - Number of Modulators Analysis on Nfkb1 : " , length( unique( cindy_analysis.env$triplets.nfkb1$Modulator ) ) )

  message("- Null Model generation for CINDy based on Mutual Information and Triplets number - Based on NFKB1")
  {
    system.time({

      mu <- mean(cindy_analysis.env$triplets.nfkb1$CMI)
      s <- sd(cindy_analysis.env$triplets.nfkb1$CMI)
      message("- Number of total rows : " , nrow(cindy_analysis.env$triplets.nfkb1) )
      message("- CMI mean : " , mu )
      message("- CMI sd : " , s )
      
      # x <- cindy_analysis.env$triplets.nfkb1 %>% group_by(Modulator) %>% summarize(n_triplets=length(Target))
      x <- cindy_analysis.env$triplets.nfkb1 %>% group_by(Modulator) %>% summarize(n_triplets=length(Target[CMI > mu]))
      y <- cindy_analysis.env$triplets.nfkb1 %>% group_by(Modulator) %>% summarize(CMI_avg_2=mean(CMI))
      
      .param <- x$n_triplets[x$n_triplets!=0]
      # .param <- x$n_triplets
      x[ x$Modulator == "Tnfrsf21" , ]
      tnf.value <- x[ x$Modulator == "Tnfrsf21" , ]$n_triplets
      
      library(fitdistrplus)
      fit.gamma <- fitdist( .param , distr = "gamma", method = "mle" )
      summary(fit.gamma)
      plot(fit.gamma)
      pgamma( tnf.value , shape = fit.gamma$estimate["shape"] , rate = fit.gamma$estimate["rate"] , lower.tail = FALSE )
      fit.pois <- fitdist( .param , distr = "pois", method = "mle" )
      summary(fit.pois)
      plot(fit.pois)
      ppois( tnf.value , lambda = fit.pois$estimate["lambda"] , lower.tail = FALSE )
      
      
      mean( x$n_triplets )
      sd( x$n_triplets )
      
      .modulators <- unique( cindy_analysis.env$triplets.nfkb1$Modulator )
      cindy_analysis.env$triplets.nfkb1[ Modulator %in% .modulators , n_triplets ]
      
      s <- sd(unique(cindy_analysis.env$triplets.nfkb1$CMI_avg))
      n <- length(unique(cindy_analysis.env$triplets.nfkb1$CMI_avg))
      message("- Number of total rows : " , n )
      message("- CMI mean : " , mu )
      message("- CMI sd : " , s )

      getCindyPvalueFromTripletsNumber <- function( triplets ) {
        .pvalue <- pgamma( triplets , shape = fit.gamma$estimate["shape"] , rate = fit.gamma$estimate["rate"] , lower.tail = FALSE )
        return(.pvalue)
      }      

    })
  }
  
  ## 
  # Function to Fix CINDy p-values bug
  # ----------------------------------
  generateCindyPvalue <- function( triplets_values , p_values ) {
    
    index.to.keep <- !( p_values == 0 )
    y <- log(p_values[ index.to.keep ])
    x <- triplets_values[ index.to.keep ]
    exponential.model <- lm(y~x)
    new_values <- triplets_values
    predicted_values <- exp( predict.lm( object = exponential.model , newdata = data.frame(x=new_values) ) )

    return(predicted_values)
  }  

  ## 
  # CINDy empirical p-value
  # -----------------------
  getCINDyEmpiricalPvalue <- function( n_triplets )
  {
    empPvals( n_triplets , cindy_analysis.env$triplets.nfkb1$n_triplets )
  }
    