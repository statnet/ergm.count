#  File R/InitWtErgmProposal.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

#' @templateVar name Disc
#' @aliases InitWtErgmProposal.Disc
#' @title Sampling for some discrete-reference ERGMs
#' @description This proposal implements [Poisson-ergmReference],
#'   [Geometric-ergmReference], [Binomial-ergmReference], and
#'   [DiscUnif-ergmReference] with arbitrary dyad level constraints.
#'
#' @concept discrete
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept valued
#'
#' @template ergmProposal-general
#'
#' @seealso [DiscTNT-ergmProposal]
NULL

InitWtErgmProposal.Disc <- function(arguments, nw, response) {
  params <- with(arguments$reference$arguments,
                 switch(arguments$reference$name,
                        Poisson = list(0L, c()),
                        Geometric = list(1L, c()),
                        Binomial = list(2L, c(trials)),
                        DiscUnif = list(3L, c(a, b))))

  list(name = "Disc", iinputs = params[[1]], inputs=as.double(params[[2]]), dyadgen = ergm_dyadgen_select(arguments, nw))
}


#' @templateVar name DiscTNT
#' @aliases InitWtErgmProposal.DiscTNT
#' @title TNT sampling for some discrete-reference ERGMs
#' @description This proposal implements [Poisson-ergmReference],
#'   [Geometric-ergmReference], [Binomial-ergmReference], and
#'   [DiscUnif-ergmReference] when the range of values includes 0,
#'   falling back to [Disc-ergmProposal] otherwise, all with arbitrary
#'   dyad-level constraints.
#'
#' @concept discrete
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept valued
#'
#' @template ergmProposal-general
#'
#' @seealso [TNT-ergmProposal]
NULL

InitWtErgmProposal.DiscTNT <- function(arguments, nw, response) {
  if(! "p0" %in% names(arguments)){
    arguments$p0 <- 0.2 
  }else if(!is.numeric(arguments$p0) || length(arguments$p0)!=1 || arguments$p0<0 || arguments$p0 >=1) {
    stop(paste("Invalid jump-to-0 argument to the TNT Discrete proposal: must be either omited or a single number in the interval [0,1)."))
  }

  params <- with(arguments$reference$arguments,
                 switch(arguments$reference$name,
                        Poisson = list(0L, c()),
                        Geometric = list(1L, c()),
                        Binomial = list(2L, c(trials)),
                        DiscUnif = list(3L, c(a, b))))

  # Fall back to Disc if DiscUnif range does not contain 0.
  if(params[[1]] == 3L && (params[[2]][1]>0 || params[[2]][2]<0)) InitWtErgmProposal.Disc(arguments, nw, response)
  else list(name = "DiscTNT", iinputs = params[[1]], inputs=as.double(c(arguments$p0, params[[2]])), dyadgen = ergm_dyadgen_select(arguments, nw))
}

#' @templateVar name ZIPoisson
#' @aliases InitWtErgmProposal.ZIPoisson
#' @title TODO
#' @description TODO
#'
#' @concept discrete
#' @concept nonnegative
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept valued
#'
#' @template ergmProposal-general
NULL

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
