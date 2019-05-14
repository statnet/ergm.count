#  File R/zzz.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.count",c("statnet"),FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage(paste(c(strwrap(paste("NOTE: The form of the term ",sQuote("CMP")," has been changed in version 3.2 of ",sQuote("ergm.count"),". See the news or help('CMP') for more information.",sep="")),""),collapse="\n"))
  }
}
  
.onLoad <- function(lib, pkg){    
  .RegisterProposals()
}

.RegisterProposals <- function(){
  ergm_proposal_table("c", "Poisson", "",  0, "random", "Poisson")
  ergm_proposal_table("c", "Poisson", "",  1, "TNT", "PoissonTNT")
  ergm_proposal_table("c", "Poisson", "",  0, "0inflated", "ZIPoisson")
  ergm_proposal_table("c", "Poisson", "observed",  0, "random", "PoissonNonObserved")

  ergm_proposal_table("c", "Geometric", "",  0, "random", "Geometric")
  ergm_proposal_table("c", "Geometric", "observed",  0, "random", "GeometricNonObserved")

  ergm_proposal_table("c", "Binomial", "",  0, "random", "Binomial")
  ergm_proposal_table("c", "Binomial", "observed",  0, "random", "BinomialNonObserved")
}
