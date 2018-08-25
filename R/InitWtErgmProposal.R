#  File R/InitWtErgmProposal.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
#######################################################################
InitWtErgmProposal.Poisson <- function(arguments, nw, response) {
  proposal <- list(name = "Poisson", inputs=NULL)
  proposal
}

InitWtErgmProposal.PoissonTNT <- function(arguments, nw, response) {
  if(! "p0" %in% names(arguments)){
    arguments$p0 <- 0.2 
  }else if(!is.numeric(arguments$p0) || length(arguments$p0)!=1 || arguments$p0<0 || arguments$p0 >=1) {
    stop(paste("Invalid jump-to-0 argument to the TNT Poisson proposal: must be either omited or a single number in the interval [0,1)."))
  }
  proposal <- list(name = "PoissonTNT", inputs=as.double(arguments$p0))
  proposal
}

InitWtErgmProposal.ZIPoisson <- function(arguments, nw, response) {
  if(! "p0" %in% names(arguments)){
    arguments$p0 <- max((1-sum(nw %e% response > 0)/network.dyadcount(nw)) -
                        exp(-sum(nw %e% response)/network.dyadcount(nw)),
                        0)
    if(arguments$p0==0) warning("The data do not appear to zero-inflated and are likely to be zero-deflated.")
    cat("Using adaptive jump-to-0 probability of ",arguments$p0,".\n")
  }else if(!is.numeric(arguments$p0) || length(arguments$p0)!=1 || arguments$p0<0 || arguments$p0 >=1) {
    stop(paste("Invalid jump-to-0 argument to the Zero-Inflated Poisson proposal: must be either omited or a single number in the interval [0,1)."))
  }
  proposal <- list(name = "ZIPoisson", inputs=as.double(arguments$p0))
  proposal
}

InitWtErgmProposal.PoissonNonObserved <- function(arguments, nw, response) {
  proposal <- list(name = "PoissonNonObserved", inputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}

InitWtErgmProposal.Binomial <- function(arguments, nw, response) {
  proposal <- list(name = "Binomial", inputs=arguments$reference$arguments$trials)
  proposal
}

InitWtErgmProposal.BinomialNonObserved <- function(arguments, nw, response) {
  proposal <- list(name = "BinomialNonObserved", inputs=c(arguments$reference$trials,to_ergm_Cdouble(is.na(nw))))
  proposal
}

InitWtErgmProposal.Geometric <- function(arguments, nw, response) {
  proposal <- list(name = "Geometric", inputs=NULL)
  proposal
}

InitWtErgmProposal.GeometricNonObserved <- function(arguments, nw, response) {
  proposal <- list(name = "GeometricNonObserved", inputs=to_ergm_Cdouble(is.na(nw)))
  proposal
}
