#  File R/ergm.count-package.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################
#' Fit, Simulate and Diagnose Exponential-Family Models for Networks with Count
#' Edges
#'
#' \code{\link[=ergm.count-package]{ergm.count}} is a set of extensions to
#' package \code{\link[=ergm-package]{ergm}} to fit and simulate from
#' exponential-family random graph models for networks whose edge weights are
#' counts. For a list of functions type \code{help(package='ergm')} and
#' \code{help(package='ergm.count')}
#'
#' Mainly, it implements Poisson, binomial, geometric, and discrete uniform
#' dyadwise reference measures for valued ERGMs
#' (documented here in [`ergmReference`]), and provides some count-specific
#' change statistics (documented in [`ergmTerm`]).
#'
#' For a complete list of the functions, use \code{library(help="ergm")} and
#' \code{library(help="ergm.count")} or read the rest of the manual.
#'
#' When publishing results obtained using this package, please cite the
#' original authors as described in \code{citation(package="ergm.count")}.
#'
#' All programs derived from this package must cite it.
#'
#' This package contains functions specific to using \code{\link{ergm}} to
#' model networks whose dyad values are counts. Examples include counts of
#' conversations, messages, and other interactions.
#'
#' In particular, this package implements the Poisson, geometric, binomial, and
#' discrete uniform reference measures (documented in
#' [`ergmReference`] for use by \code{\link{ergm}} and
#' \code{\link{simulate.ergm}}) to fit models from this family, as well as
#' statistics specific to modeling counts, such as the [`CMP`][CMP-ergmTerm] for
#' the Conway-Maxwell-Poisson Distribution.
#'
#' For detailed information on how to download and install the software, go to
#' the Statnet project website: \url{https://statnet.org}. A tutorial, support
#' newsgroup, references and links to further resources are provided there.
#'
#' @section Known issues:
#'
#' ## Parameter space constraints
#'
#' Poisson- and geometric-reference ERGMs have an unbouded sample space. This
#' means that the parameter space may be constrained in complex ways that
#' depend on the terms used in the model. At this time [`ergm`] has
#' no way to detect when a parameter configuration had strayed outside of the
#' parameter space, but it may be noticeable on a runtime trace plot (activated
#' via `MCMC.runtime.traceplot` control parameter), when the simulated
#' values keep climbing upwards. (See Krivitsky (2012) for a further
#' discussion.)
#'
#' A possible remedy if this appears to occur is to try lowering the control
#' parameter `MCMLE.steplength`.
#'
#' @name ergm.count-package
#' @docType package
#' @author Pavel N. Krivitsky \email{pavel@@statnet.org}
#' @seealso [`ergmTerm`], [`ergmReference`]
#' @references Handcock MS, Hunter DR, Butts CT, Goodreau SG, Krivitsky PN and
#' Morris M (2012). _Fit, Simulate and Diagnose Exponential-Family Models for
#' Networks_. Version 3.1.  Project home page at <URL: https://www.statnet.org>,
#' <URL: CRAN.R-project.org/package=ergm>.
#'
#' Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' Shmueli G, Minka TP, Kadane JB, Borle S, and Boatwright P (2005). A
#' Useful Distribution for Fitting Discrete Data: Revival of the
#' Conway--Maxwell--Poisson Distribution. *Journal of the Royal
#' Statistical Society: Series C*, 54(1): 127-142.
#'
#' @keywords package models
NULL
