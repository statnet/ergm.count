InitWtMHP.Poisson <- function(arguments, nw, response) {
  MHproposal <- list(name = "Poisson", inputs=NULL)
  MHproposal
}

InitWtMHP.PoissonTNT <- function(arguments, nw, response) {
  if(! "p0" %in% names(arguments)){
    arguments$p0 <- 0.2 
  }else if(!is.numeric(arguments$p0) || length(arguments$p0)!=1 || arguments$p0<0 || arguments$p0 >=1) {
    stop(paste("Invalid jump-to-0 argument to the TNT Poisson proposal: must be either omited or a single number in the interval [0,1)."))
  }
  MHproposal <- list(name = "PoissonTNT", inputs=as.double(arguments$p0))
  MHproposal
}

InitWtMHP.ZIPoisson <- function(arguments, nw, response) {
  if(! "p0" %in% names(arguments)){
    arguments$p0 <- max((1-sum(nw %e% response > 0)/network.dyadcount(nw)) -
                        exp(-sum(nw %e% response)/network.dyadcount(nw)),
                        0)
    if(arguments$p0==0) warning("The data do not appear to zero-inflated and are likely to be zero-deflated.")
    cat("Using adaptive jump-to-0 probability of ",arguments$p0,".\n")
  }else if(!is.numeric(arguments$p0) || length(arguments$p0)!=1 || arguments$p0<0 || arguments$p0 >=1) {
    stop(paste("Invalid jump-to-0 argument to the Zero-Inflated Poisson proposal: must be either omited or a single number in the interval [0,1)."))
  }
  MHproposal <- list(name = "ZIPoisson", inputs=as.double(arguments$p0))
  MHproposal
}

InitWtMHP.PoissonNonObserved <- function(arguments, nw, response) {
  MHproposal <- list(name = "PoissonNonObserved", inputs=ergm.Cprepare.miss(nw))
  MHproposal
}

InitWtMHP.Binomial <- function(arguments, nw, response) {
  MHproposal <- list(name = "Binomial", inputs=arguments$reference$trials)
  MHproposal
}

InitWtMHP.BinomialNonObserved <- function(arguments, nw, response) {
  MHproposal <- list(name = "BinomialNonObserved", inputs=c(arguments$reference$trials,ergm.Cprepare.miss(nw)))
  MHproposal
}

InitWtMHP.Geometric <- function(arguments, nw, response) {
  MHproposal <- list(name = "Geometric", inputs=NULL)
  MHproposal
}

InitWtMHP.GeometricNonObserved <- function(arguments, nw, response) {
  MHproposal <- list(name = "GeometricNonObserved", inputs=ergm.Cprepare.miss(nw))
  MHproposal
}
