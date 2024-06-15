#  File R/zzz.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################
#' @import statnet.common network ergm
#' @useDynLib ergm.count
#' @importFrom Rdpack reprompt
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
  ergm_proposal_table("c", c("Poisson","Binomial","Geometric","DiscUnif"), "|.dyads",  0, "random", "Disc")
  ergm_proposal_table("c", c("Poisson","Binomial","Geometric","DiscUnif"), "&sparse|.dyads",  1, "TNT", "DiscTNT")
  ergm_proposal_table("c", "Poisson", "&sparse",  0, "0inflated", "ZIPoisson")
}
