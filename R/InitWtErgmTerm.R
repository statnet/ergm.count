#  File R/InitWtErgmTerm.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

#' @templateVar name CMP
#' @title Conway-Maxwell-Poisson Distribution
#' @description This term adds one statistic to the model, of the form
#'   \eqn{\sum_{i,j}\log(y_{i,j}!)} . This turns a Poisson- or a
#'   geometric-reference ERGM into a Conway-Maxwell-Poisson-reference
#'   ERGM, allowing it to represent a broad range of disperson
#'   values. In particular, combined with a Poisson-reference ERGM, a
#'   negative coefficient on this term induces underdispersion and a
#'   positive coefficient induces overdispersion. (This behavior is
#'   different from 3.1.1, when the negation of this value was used.)
#'
#' @details Note that its current implementation may not perform well if the
#'   data are overdispersed relative to geometric.
#'
#' @usage
#' # valued: CMP
#'
#' @concept nonnegative
#' @concept directed
#' @concept undirected
#'
#' @template ergmTerm-general
InitWtErgmTerm.CMP<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="sumlogfactorial",
       coef.names="CMP",
       inputs=NULL,
       dependence=FALSE,
       minval=0)
}

#' @templateVar name CMB
#' @title Conway-Maxwell-Binomial Distribution
#' @description If `couple==TRUE` , this
#'   term adds one statistic to the model, of the form
#'   \eqn{\sum_{i,j}\log(y_{i,j}!) + \log(t-y_{i,j}!)} . This turns a Binomial- or a
#'   discrete-uniform-reference ERGM into a Conway-Maxwell-Binomial-reference
#'   ERGM, allowing it to represent a broad range of disperson
#'   values. In particular, combined with a Binomial-reference ERGM, a
#'   negative coefficient on this term induces underdispersion and a
#'   positive coefficient induces overdispersion.
#'
#'   If `coupled==FALSE` the two summands above are added as their own
#'   statistic (each with its own free parameter).
#'
#' @usage
#' # valued: CMB(trials, coupled = TRUE)
#' @param trails model parameter
#' @param coupled logical
#'
#' @concept nonnegative
#' @concept directed
#' @concept undirected
#'
#' @template ergmTerm-general
InitWtErgmTerm.CMB<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("trials", "coupled"),
                      vartypes = c("numeric", "logical"),
                      defaultvalues = list(NULL, TRUE),
                      required = c(TRUE, FALSE) )
  list(name=if(a$coupled) "CMB" else "CMB2",
       coef.names=if(a$coupled) "CMB" else c("CMB.y","CMB.n-y"),
       inputs=a$trials,
       dependence=FALSE,
       minval=0)
}
