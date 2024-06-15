#  File R/ergm.count-package.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################
#' Fit, Simulate and Diagnose Exponential-Family Models for Networks with Count
#' Edges
#'
#' \code{\link[=ergm.count-package]{ergm.count}} is a set of extensions to
#' package \code{\link[=ergm-package]{ergm}} to fit and simulate from
#' exponential-family random graph models for networks whose edge weights are
#' counts \insertCite{Kr12e}{ergm.count}.
#'
#' Mainly, it implements Poisson, binomial, geometric, and discrete
#' uniform dyadwise reference measures for valued ERGMs (documented
#' here in [`ergmReference`]), and provides some count-specific change
#' statistics (documented in [`ergmTerm`])
#' \insertCite{Kr12e,KrHu23e}{ergm.count}, including
#' [`CMP`][CMP-ergmTerm] for the Conway--Maxwell--Poisson Distribution
#' \insertCite{ShMi05u}{ergm.count}.

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
#' values keep climbing upwards. (See \insertCite{Kr12e;textual}{ergm.count} for a further
#' discussion.)
#'
#' A possible remedy if this appears to occur is to try lowering the
#' [control.ergm()] parameter `MCMLE.steplength`.
#'
#' @name ergm.count-package
#' @seealso [`ergmTerm`], [`ergmReference`]
#' @references \insertAllCited{}
"_PACKAGE"
