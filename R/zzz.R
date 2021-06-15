#  File R/zzz.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################
#' @import statnet.common network ergm
#' @useDynLib ergm.count
.onAttach <- function(libname, pkgname){
  sm <- statnetStartupMessage("ergm.count",c("statnet"),FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}
  
.onLoad <- function(libname, pkgname){
  .RegisterProposals()
}

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Poisson", "",  0, "random", "Poisson")
  ergm_proposal_table("c", "Poisson", "&sparse",  1, "TNT", "PoissonTNT")
  ergm_proposal_table("c", "Poisson", "&sparse",  0, "0inflated", "ZIPoisson")
  ergm_proposal_table("c", "Poisson", "observed",  0, "random", "PoissonNonObserved")

  ergm_proposal_table("c", "Geometric", "",  0, "random", "Geometric")
  ergm_proposal_table("c", "Geometric", "observed",  0, "random", "GeometricNonObserved")

  ergm_proposal_table("c", "Binomial", "",  0, "random", "Binomial")
  ergm_proposal_table("c", "Binomial", "observed",  0, "random", "BinomialNonObserved")
}
