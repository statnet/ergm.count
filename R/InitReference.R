#  File R/InitErgmReference.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2016 Statnet Commons
#######################################################################
InitErgmReference.Poisson <- function(lhs.nw, ...){
  list(name="Poisson")
}

InitErgmReference.Binomial <- function(lhs.nw, trials, ...){
  
  list(name="Binomial", trials=trials)
}

InitErgmReference.Geometric <- function(lhs.nw, ...){
  
  list(name="Geometric")
}
