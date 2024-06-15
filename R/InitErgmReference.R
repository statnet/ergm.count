#  File R/InitErgmReference.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################

#' @templateVar name Poisson
#' @title Poisson-reference ERGM
#' @description Specifies each
#'   dyad's baseline distribution to be Poisson with mean 1:
#'   \eqn{h(y)=\prod_{i,j} 1/y_{i,j}!} , with the support of
#'   \eqn{y_{i,j}} being natural numbers (and \eqn{0} ). Using
#'   [valued ERGM terms][ergmTerm] that are
#'   "generalized" from their binary counterparts, with form
#'   `"sum"` (see previous link for the list) produces Poisson
#'   regression. Using [`CMP`][CMP-ergmTerm] induces a
#'   Conway-Maxwell-Poisson distribution that is Poisson when its
#'   coefficient is \eqn{0} and geometric when its coefficient is
#'   \eqn{1} .
#'
#'  @details Three proposal functions are currently implemented, two of them
#'   designed to improve mixing for sparse networks. They can can be
#'   selected via the `MCMC.prop.weights=` control parameter. The
#'   sparse proposals work by proposing a jump to 0. Both of them take
#'   an optional proposal argument `p0` (i.e.,
#'   `MCMC.prop.args=list(p0=...)` ) specifying the probability of
#'   such a jump. However, the way in which they implement it are
#'   different:
#'
#'   - `"random"`: Select a dyad (i,j) at random, and draw the
#'   proposal \eqn{y_{i,j}^\star \sim \mathrm{Poisson}_{\ne
#'   y_{i,j}}(y_{i,j}+0.5)} (a Poisson distribution with mean
#'   slightly higher than the current value and conditional on not
#'   proposing the current value).
#'
#'   - `"0inflated"`: As `"random"` but, with
#'   probability `p0` , propose a jump to 0 instead of a
#'   Poisson jump (if not already at 0). If `p0` is not given,
#'   defaults to the "surplus" of 0s in the observed network,
#'   relative to Poisson.
#'
#'   - `"TNT"`: (the default) As `"0inflated"` but
#'   instead of selecting a dyad at random, select a tie with
#'   probability `p0` , and a random dyad otherwise, as with
#'   the binary TNT. Currently, `p0` defaults to 0.2.
#'
#' @usage
#' # Poisson
#' @concept discrete
#' @concept nonnegative
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept valued
#'
#' @template ergmReference-general
InitErgmReference.Poisson <- function(lhs.nw, ...){
  list(name="Poisson", init_methods = c("CD","zeros"))
}

#' @templateVar name Binomial
#' @title Binomial-reference ERGM
#' @description Specifies each dyad's baseline distribution to be binomial with
#'   `trials` trials and success probability of \eqn{0.5} :
#'   \eqn{h(y)=\prod_{i,j}{{\mathrm{trials}}\choose{y_{i,j}}}} . Using
#'   [valued ERGM terms][ergmTerm] that are
#'   "generalized" from their binary counterparts, with form
#'   `"sum"` (see previous link for the list) produces logistic
#'   regression.
#'
#' @usage
#' # Binomial(trials)
#' @param trails model parameter
#' @concept discrete
#' @concept finite
#' @concept nonnegative
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept valued
#'
#' @template ergmReference-general
InitErgmReference.Binomial <- function(lhs.nw, trials, ...){
  list(name="Binomial", arguments=list(trials=trials), init_methods = c("CD","zeros"))
}

#' @templateVar name Geometric
#' @title Geometric-reference ERGM
#' @description Specifies
#'   each dyad's baseline distribution to be uniform on the natural
#'   numbers (and \eqn{0} ): \eqn{h(y)=1} . In itself, this
#'   "distribution" is improper, but in the presence of
#'   [`sum`][sum-ergmTerm] , a geometric
#'   distribution is induced. Using [`CMP`][CMP-ergmTerm] (in addition to
#'   [`sum`][sum-ergmTerm] ) induces a
#'   Conway-Maxwell-Poisson distribution that is geometric when its
#'   coefficient is \eqn{0} and Poisson when its coefficient is
#'   \eqn{-1} .
#'
#' @usage
#' # Geometric
#'
#' @concept discrete
#' @concept nonnegative
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept valued
#'
#' @template ergmReference-general
InitErgmReference.Geometric <- function(lhs.nw, ...){
  list(name="Geometric", init_methods = c("CD","zeros"))
}
