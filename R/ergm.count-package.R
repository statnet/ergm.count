#  File R/ergm.count-package.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
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
#' (\link[=ergm-references]{documented here}), and provides some count-specific
#' change statistics (\link[=ergm-terms]{documented here}).
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
#' \code{\link{ergm-references}}) for use by \code{\link{ergm}} and
#' \code{\link{simulate.ergm}} to fit models from this family, as well as
#' statistics specific to modeling counts, such as the \code{\link{CMP}} for
#' the Conway-Maxwell-Poisson Distribution.
#'
#' For detailed information on how to download and install the software, go to
#' the Statnet project website: \url{https://statnet.org}. A tutorial, support
#' newsgroup, references and links to further resources are provided there.
#'
#' @name ergm.count-package
#' @docType package
#' @author Pavel N. Krivitsky \email{pavel@@statnet.org}
#' @seealso \code{\link{ergm-terms}}, \code{\link{ergm-references}}
#' @references Handcock MS, Hunter DR, Butts CT, Goodreau SG, Krivitsky PN and
#' Morris M (2012). _Fit, Simulate and Diagnose Exponential-Family Models for
#' Networks_. Version 3.1.  Project home page at <URL: https://www.statnet.org>,
#' <URL: CRAN.R-project.org/package=ergm>.
#'
#' Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#' @keywords package models
NULL
