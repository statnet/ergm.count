#  File R/InitErgmReference.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
#######################################################################
InitErgmReference.Poisson <- function(lhs.nw, ...){
  list(name="Poisson", init_methods = c("CD","zeros"))
}

InitErgmReference.Binomial <- function(lhs.nw, trials, ...){
  list(name="Binomial", arguments=list(trials=trials), init_methods = c("CD","zeros"))
}

InitErgmReference.Geometric <- function(lhs.nw, ...){
  list(name="Geometric", init_methods = c("CD","zeros"))
}
