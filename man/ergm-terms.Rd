\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{CMP}
\alias{InitWtErgmTerm.CMP}
\docType{methods}
\title{Terms used in Exponential Family Random Graph Models Specific to Counts}
\description{
  This page describes the possible terms (and hence network statistics)
  specific for modeling count data.

  See the \code{\link[ergm]{ergm-terms}} documentation in the
  \code{ergm} package for a complete description of what ERGM terms are
  and how they are used.
}

\section{Possible terms to represent network statistics}{

  The terms currently available are:
  \describe{
    
    \item{\code{CMP}}{\emph{Conway-Maxwell-Poisson Distribution:}
      This term adds one statistic to the model, of the form
      \eqn{\sum_{i,j}-\log(y_{i,j}!)}. This turns a Poisson- or a
      geometric-reference ERGM into a Conway-Maxwell-Poisson-reference
      ERGM, allowing it to represent a broad range of disperson values.

      Note that its current implementation may not perform well if the
      data are overdispersed relative to geometric.
    }
  }
}
\references{
Handcock MS, Hunter DR, Butts CT, Goodreau SG, Krivitsky PN and Morris M
(2012). _Fit, Simulate and Diagnose Exponential-Family Models for
Networks_. Version 3.1.  Project home
page at <URL: http://www.statnet.org>, <URL:
CRAN.R-project.org/package=ergm>.
  
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}

Shmueli G, Minka TP, Kadane JB, Borle S, and Boatwright P (2005). A Useful Distribution for Fitting Discrete Data: Revival of the Conway--Maxwell--Poisson Distribution. \emph{Journal of the Royal Statistical Society: Series C}, 54(1): 127-142.

}
\seealso{ergm, network, \%v\%, \%n\%, sna,
         summary.ergm, print.ergm}
\keyword{models}
