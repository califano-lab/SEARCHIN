##
# Utility Functions for VIPER
# ----------------------------------------------
# @author: Mariano Alvarez
# @copyright: 2019
# ----------------------------------------------

#' Variance stabilization transformation for RNAseq data
#' This function stabilizes the variance, transform the data and add shot noise to the data
#' @param x CountDataSet or matrix containing the raw counts, with genes in rows and samples in columns
#' @param method Character string indicating the method for estimatinf the dispersion (see DESeq::estimateDispersions)
#' @param fitType Character string indicating the type of fit for the dispersion (see DESeq::estimateDispersions)
#' @return Expression matrix
#' @export

DEtransform <- function(x, method=c("blind", "pooled", "pooled-CR", "per-condition"), fitType=c("parametric", "local"), noise=TRUE) {
  method <- match.arg(method)
  fitType <- match.arg(fitType)
  cnames <- NULL
  if (class(x) != "CountDataSet") {
    if (class(x) != "matrix") stop("x must be a CountDataSet or integer matrix object", call.=F)
    if (length(which(duplicated(colnames(x))))>0) {
      cnames <- colnames(x)
      colnames(x) <- 1:ncol(x)
    } 
    x <- newCountDataSet(x, factor(colnames(x)))
  }
  x <- estimateSizeFactors(x)
  x <- estimateDispersions(x, method=method, fitType=fitType)
  x <- getVarianceStabilizedData(x)
  tmp <- x
  if (noise) {
    tmp <- apply(x, 2, function(x) {
      x <- sort(unique(x))
      x <- cbind(x[1:(length(x)-1)], x[2:length(x)])
      x <- cbind(x[, 1], sqrt(rvar(x)))
      return(x)
    })
    tmp <- cbind(unlist(lapply(tmp, function(x) x[, 1]), use.names=F), unlist(lapply(tmp, function(x) x[, 2]), use.names=F))
    tmp1 <- smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
    tmp[tmp[, 1]>tmp1$x[which.min(tmp1$y)], 2] <- 0
    tmp1 <- smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
    tmp <- x+rnorm(length(x))*predict(tmp1, x)$y
  }
  if (!is.null(cnames)) colnames(tmp) <- cnames
  return(tmp)
}

#' Variance of rows for arrays with NA values
#' This function computes the variance by rows ignoring NA values
#' @param x Numeric matrix
#' @return Numeric vector
#' @export


rvar <- function(x) {
  ave <- as.vector(rmean(x))
  pos <- which(is.na(x))
  largo <- rlength(x)
  x[pos] <- rep(ave,ncol(x))[pos]
  (x-ave)^2 %*% rep(1,ncol(x))/(largo-1)
}

#' Mean of rows for arrays with NA values
#' This function compute the mean by rows ignoring NA values
#' @param x Numeric matrix
#' @return Numeric vector
#' @export

rmean <- function(x) {
  largo <- rlength(x)
  x[is.na(x)] <- 0
  res <- x %*% rep(1,ncol(x)) / largo
  names(res) <- rownames(x)
  res
}


#' Length of rows for arrays with NA values
#' This function report the length of rows of a matrix ignoring the NA values
#' @param x Matrix
#' @return Integer indicating the number of non-NA rows
#' @export

rlength <- function(x) {
  r <- x/x
  r[x==0] <- 1
  r[!is.finite(r)] <- 0
  r %*% rep(1,ncol(r))
}

# filter for expressed genes in RNAseq data
#' RNAseq expression filter
#' This function remove the non-expressed genes based on fitting a mixture of two gaussians to the distribution of the normalized (DEtransform) RNAseq signal
#' @param dset DEtransform normalized expression matrix
#' @param dist Integer indicating the number of distributions to fit, usually 2 or 3
#' @param output Character string indicating the output type, either likelihood or dset
#' @param threshold Number indicating the threhold for the relative likelihood of expression
#' @param method Character string indicating whether the analysis should be performed sample-by-sample (sample) or globally (global)
#' @param maxsize Integer indicating the maximum size for the dataset to fit the Gaussian mixture model if global method is selected
#' @param cores Integer indicating the number of cores to use
#' @return Filtered expression matrix
#' @export

rnaseqExpressionFilter <- function(dset, dist=2, output=c("likelihood", "dset"), threshold=.8, method=c("sample", "global"), maxsize=1e6, cores=1) {
  output <- match.arg(output)
  method <- match.arg(method)
  switch(method,
         sample={
           rl <- mclapply(split(t(dset), 1:ncol(dset)), function(x) {
             tmp <- get2peaks(x)
             pos <- seq(tmp[1], tmp[2], length=dist)[-c(1, dist)]
             fit <- list(mu = c(tmp, pos), sigma = c(1, 2, rep(.5, length(pos))), lambda = rep(.5, dist), all.loglik = rep(0, 1001))
             count <- 0
             x1 <- density(x, na.rm=TRUE)
             x1 <- rep(x1$x, round(x1$y*1e5/sum(x1$y)))
             while (length(fit$all.loglik) > 1000 & count < 3) {
               fit <- normalmixEM(x1, mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda, maxit = 1000, verb = FALSE)
               count <- count + 1
             }
             a <- pnorm(x, mean=fit$mu[2], sd=fit$sigma[2], lower.tail=TRUE)
             b <- pnorm(x, mean=fit$mu[1], sd=fit$sigma[1], lower.tail=FALSE)
             a/(a+b)
           }, mc.cores=cores)
           rl <- do.call(cbind, rl)
           colnames(rl) <- colnames(dset)
           rownames(rl) <- rownames(dset)
         },
         global={
           tmp <- get2peaks(dset)
           pos <- seq(tmp[1], tmp[2], length=dist)[-c(1, dist)]
           fit <- list(mu = c(tmp, pos), sigma = c(1, 2, rep(.5, length(pos))), lambda = rep(.5, dist), all.loglik = rep(0, 1001))
           x1 <- density(dset, na.rm=TRUE)
           x1 <- rep(x1$x, round(x1$y*1e5/sum(x1$y)))
           count <- 0
           while (length(fit$all.loglik) > 1000 & count < 3) {
             fit <- normalmixEM(x1, mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda, maxit = 1000, verb = FALSE)
             count <- count + 1
           }
           a <- pnorm(dset, mean=fit$mu[2], sd=fit$sigma[2], lower.tail=TRUE)
           b <- pnorm(dset, mean=fit$mu[1], sd=fit$sigma[1], lower.tail=FALSE)
           rl <- a/(a+b)
         })
  if (output=="likelihood") return(rl)
  return(filterRowMatrix(dset, rowSums(rl<threshold)==0))
}

get2peaks <- function(x) {
  den <- density(x, adj=1.5, na.rm=TRUE)
  posi <- which(c(FALSE, diff(diff(den$y)>0)!=0, FALSE))
  posi <- posi[order(den$y[posi], decreasing=TRUE)[1:2]]
  sort(den$x[posi])
}

#MAP entrezID to gene
entrez2gene <- function(id) {
  if (!exists("desc", envir=as.environment(".GlobalEnv"))) data(desc, envir=as.environment(".GlobalEnv"))
  desc <- get("desc", env=as.environment(".GlobalEnv"))
  res <- desc[match(id, desc[, 2]), 3]
  names(res) <- id
  res
}
