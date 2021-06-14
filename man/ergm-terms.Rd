%  File man/ergm-terms.Rd in package ergm.count, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2008-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{CMP}
\alias{CMB}
\alias{InitWtErgmTerm.CMP}
\alias{InitWtErgmTerm.CMB}
\docType{methods}
\title{Terms used in Exponential Family Random Graph Models Specific to Counts}
\description{
  This page describes the possible terms (and hence network statistics)
  included in the \code{\link[=ergm.count-package]{ergm.count}} package.

  See the \code{\link[ergm]{ergm-terms}} documentation in the
  \code{ergm} package for a complete description of what ERGM terms are
  and how they are used.
}

\section{Terms to represent network statistics included in the
  \code{\link[=ergm.count-package]{ergm.count}} pacakge}{

  All terms listed are valued.
  
  \describe{
    
    \item{\code{CMP}}{\emph{Conway-Maxwell-Poisson Distribution:} This
      term adds one statistic to the model, of the form
      \eqn{\sum_{i,j}\log(y_{i,j}!)}. This turns a Poisson- or a
      geometric-reference ERGM into a Conway-Maxwell-Poisson-reference
      ERGM, allowing it to represent a broad range of disperson
      values. In particular, combined with a Poisson-reference ERGM, a
      negative coefficient on this term induces underdispersion and a
      positive coefficient induces overdispersion. (This behavior is
      different from 3.1.1, when the negation of this value was used.)

      Note that its current implementation may not perform well if the
      data are overdispersed relative to geometric.
    }

    \item{\code{CMB(trials, coupled = TRUE)}}{
      \emph{Conway-Maxwell-Binomial Distribution:} If
      \code{couple==TRUE}, this
      term adds one statistic to the model, of the form
      \eqn{\sum_{i,j}\log(y_{i,j}!) + \log(t-y_{i,j}!)}. This turns a Binomial- or a
      discrete-uniform-reference ERGM into a Conway-Maxwell-Binomial-reference
      ERGM, allowing it to represent a broad range of disperson
      values. In particular, combined with a Binomial-reference ERGM, a
      negative coefficient on this term induces underdispersion and a
      positive coefficient induces overdispersion.

      If \code{coupled==FALSE} the two summands above are added as their own
      statistic (each with its own free parameter).
    }
  }
}
\references{
  \itemize{
    \item Handcock M. S., Hunter D. R., Butts C. T., Goodreau S. G.,
      Krivitsky P. N. and Morris M.
(2012). _Fit, Simulate and Diagnose Exponential-Family Models for
Networks_. Version 3.1.  Project home
page at <URL: https://www.statnet.org>, <URL:
CRAN.R-project.org/package=ergm>.
  
\item Krivitsky P. N. (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \doi{10.1214/12-EJS696}

\item Shmueli G., Minka T. P., Kadane J. B., Borle S., and Boatwright P. (2005). A Useful Distribution for Fitting Discrete Data: Revival of the Conway--Maxwell--Poisson Distribution. \emph{Journal of the Royal Statistical Society: Series C}, 54(1): 127-142.

\item Shmueli G., Minka T. P., Kadane J. B., Borle S., and Boatwright P. (2005). A Useful Distribution for Fitting Discrete Data: Revival of the Conway--Maxwell--Poisson Distribution. \emph{Journal of the Royal Statistical Society: Series C}, 54(1): 127-142.

\item Kadane, Joseph B.  (2016) Sums of Possibly Associated Bernoulli Variables: The Conway-Maxwell-Binomial Distribution. \emph{Bayesian Analysis}, 11(2): 403--420. \doi{10.1214/15-BA955}

}
}
\seealso{\code{\link[ergm:ergm-terms]{ergm-terms}} (from the \code{\link[=ergm-package]{ergm}} package), \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\keyword{models}
