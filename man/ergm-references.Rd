\name{ergm-references}
\alias{ergm-references}
\alias{references-ergm}
\alias{ergm.references}
\alias{references.ergm}
\alias{Poisson}
\alias{Binomial}
\alias{DiscUnif}
\alias{Geometric}
\alias{InitReference.Poisson}
\alias{InitReference.Binomial}
\alias{InitReference.DiscUnif}
\alias{InitReference.Geometric}
\alias{InitWtMHP.Poisson}
\alias{InitWtMHP.ZIPoisson}
\alias{InitWtMHP.PoissonNonObserved}
\alias{InitWtMHP.Binomial}
\alias{InitWtMHP.BinomialNonObserved}
\alias{InitWtMHP.Geometric}
\alias{InitWtMHP.GeometricNonObserved}
\alias{InitWtMHP.DiscUnif}
\alias{InitWtMHP.DiscUnifNonObserved}

\docType{methods}

\title{Reference Measures for Exponential-Family Random Graph Models for
  Counts}

\description{
  This page describes the possible reference measures (baseline distributions)
  for modeling count data.

  Each of these is specified on the RHS of a one-sided formula passed as
  the \code{reference} argument to \code{\link[ergm]{ergm}}. 
  See the \code{\link[ergm]{ergm}} documentation for a complete
  description of how reference measures are specified.
}

\section{Possible reference measures to represent baseline
  distributions}{

  Reference measures currently available are:
  \describe{
    
    \item{\code{Poisson}}{\emph{Poisson-reference ERGM:} Specifies each
      dyad's baseline distribution to be Poisson with mean 1:
      \eqn{h(y)=\prod_{i,j} 1/y_{i,j}!}, with the support of
      \eqn{y_{i,j}} being natural numbers (and \eqn{0}). Using
      \code{\link[ergm:ergm-terms]{valued ERGM terms}} that are
      ``generalized'' from their binary counterparts, with form
      \code{"sum"} (see previous link for the list) produces Poisson
      regression. Using \code{\link{CMP}} induces a
      Conway-Maxwell-Poisson distribution that is Poisson when its
      coefficient is \eqn{0} and geometric when its coefficient is
      \eqn{-1}.

      A 0-inflated version of this proposal is also implemented. It can
      be invoked by specifying a \code{\link[ergm]{control.ergm}} (or
      \code{\link[ergm]{control.simulate.formula}}, etc.)  parameter
      \code{MCMC.prop.weights="0inflated"}. An optional proposal
      argument \code{p0} (i.e., \code{MCMC.prop.args=list(p0=...)}) can
      be passed as well, specifying the probability of a jump to 0. (If
      not given, defaults to the "surplus" of 0s in the observed
      network, relative to Poisson.)  }

    \item{\code{Geometric}}{\emph{Geometric-reference ERGM:} Specifies
      each dyad's baseline distribution to be uniform on the natural
      numbers (and \eqn{0}): \eqn{h(y)=1}. In itself, this
      "distribution" is improper, but in the presence of
      \code{\link[ergm:ergm-terms]{sum}}, a geometric
      distribution is induced. Using \code{\link{CMP}} (in addition to
      \code{\link[ergm:ergm-terms]{sum}}) induces a
      Conway-Maxwell-Poisson distribution that is geometric when its
      coefficient is \eqn{0} and Poisson when its coefficient is
      \eqn{1}.  }

    \item{\code{Binomial(trials)}}{\emph{Binomial-reference ERGM:}
      Specifies each dyad's baseline distribution to be binomial with
      \code{trials} trials and success probability of \eqn{0.5}:
      \eqn{h(y)=\prod_{i,j}\choose{trials}{y_{i,j}}}.  Using
      \code{\link[ergm:ergm-terms]{valued ERGM terms}} that are
      ``generalized'' from their binary counterparts, with form
      \code{"sum"} (see previous link for the list) produces logistic
      regression.  }

    \item{\code{DiscUnif(a,b)}}{\emph{Discrete-Uniform-reference ERGM:}
      Specifies each dyad's baseline distribution to be discrete uniform
      between \code{a} and \code{b} (both inclusive): \eqn{h(y)=1}, with
      the support being
      \eqn{\code{a},\code{a}+1,\ldots,\code{b}-1,\code{b}}.  }

    }
}

\references{
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}

Shmueli G, Minka TP, Kadane JB, Borle S, and Boatwright P (2005). A Useful Distribution for Fitting Discrete Data: Revival of the Conway--Maxwell--Poisson Distribution. \emph{Journal of the Royal Statistical Society: Series C}, 54(1): 127-142.
}

\seealso{ergm, network, \%v\%, \%n\%, sna,
         summary.ergm, print.ergm}
\keyword{models}