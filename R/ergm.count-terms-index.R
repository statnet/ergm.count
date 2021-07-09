#  File R/ergm-count-terms-index.R in package ergm.count, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2021 Statnet Commons
#######################################################################

#' Terms used in Exponential Family Random Graph Models Specific to Counts
#'
#' @name ergmTerm
#' @aliases ergm-terms ergm.terms terms-ergm terms.ergm
#' @docType package
#' @description This page describes the possible terms (and hence network statistics)
#' included in the [`ergm.count`][ergm.count-package] package.
#'
#' See the \code{\link[ergm:ergm-terms]{ergm-terms}} documentation in the `ergm`
#' package for a complete description of what ERGM terms are and how they are
#' used.
#'
#' @section Terms to represent network statistics included in the [`ergm.count`][ergm.count-package] package:
#' ## Term index
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmTerm"))}}
#'
#' ## Frequently-used terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::FREQUENTLY_USED_TERM_CATEGORIES, only.include='frequently-used'))}}
#'
#' ## Operator terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", categories=ergm:::OPERATOR_CATEGORIES, only.include='operator'))}}
#'
#' ## All terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm"))}}
#'
#' ## Terms by concept
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmTerm"))}}
#'
#' @seealso \code{\link[ergm:ergm-terms]{ergm-terms}} (from the
#' [`ergm`][ergm-package] package), [`ergm`], [`network`], `%v%`, `%n%`
#'
#' @references
#' - Handcock M. S., Hunter D. R., Butts C. T.,
#' Goodreau S. G., Krivitsky P. N. and Morris M. (2012). _Fit, Simulate and
#' Diagnose Exponential-Family Models for Networks_. Version 3.1.  Project home
#' page at <URL: https://www.statnet.org>, <URL:
#' CRAN.R-project.org/package=ergm>.
#'
#' - Krivitsky P. N. (2012). Exponential-Family Random Graph Models for
#' Valued Networks. *Electronic Journal of Statistics*, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#'
#' - Shmueli G., Minka T. P., Kadane J. B., Borle S., and Boatwright P.
#' (2005). A Useful Distribution for Fitting Discrete Data: Revival of the
#' Conway--Maxwell--Poisson Distribution. *Journal of the Royal
#' Statistical Society: Series C*, 54(1): 127-142.
#'
#' - Shmueli G., Minka T. P., Kadane J. B., Borle S., and Boatwright P.
#' (2005). A Useful Distribution for Fitting Discrete Data: Revival of the
#' Conway--Maxwell--Poisson Distribution. *Journal of the Royal
#' Statistical Society: Series C*, 54(1): 127-142.
#'
#' - Kadane, Joseph B.  (2016) Sums of Possibly Associated Bernoulli
#' Variables: The Conway-Maxwell-Binomial Distribution. *Bayesian
#' Analysis*, 11(2): 403--420. \doi{10.1214/15-BA955}
#'
#' @keywords models
NULL

#' Reference Measures for Exponential-Family Random Graph Models for Counts
#'
#' @name ergmReference
#' @aliases ergm-references references-ergm ergm.references references.ergm Poisson Binomial Geometric InitWtErgmProposal.Poisson InitWtErgmProposal.ZIPoisson InitWtErgmProposal.PoissonTNT InitWtErgmProposal.PoissonNonObserved InitWtErgmProposal.Binomial InitWtErgmProposal.BinomialNonObserved InitWtErgmProposal.Geometric InitWtErgmProposal.GeometricNonObserved
#' @docType package
#' @description This page describes the possible reference measures (baseline distributions)
#' for modeling count data found in the [`ergm.count`][ergm.count-package] package.
#'
#' Each of these is specified on the RHS of a one-sided formula passed as the
#' `reference` argument to [`ergm`].  See the [`ergm`] documentation for a complete description of how
#' reference measures are specified.
#'
#' @section Known issues:
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
#' @section Possible reference measures to represent baseline distributions:
#' ## Reference index
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmReference"))}}
#'
#' @seealso [`ergm`][ergm-package], [`network`], `\%v\%`, `\%n\%`, `sna`, [`summary.ergm`], [`print.ergm`]
#'
#' @references
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for
#' Valued Networks. *Electronic Journal of Statistics*, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#'
#' - Shmueli G, Minka TP, Kadane JB, Borle S, and Boatwright P (2005). A Useful
#' Distribution for Fitting Discrete Data: Revival of the
#' Conway--Maxwell--Poisson Distribution. *Journal of the Royal
#' Statistical Society: Series C*, 54(1): 127-142.
#'
#' @keywords models
NULL
