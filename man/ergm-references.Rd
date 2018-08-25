%  File man/ergm-references.Rd in package ergm.count, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2008-2018 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-references}
\alias{ergm-references}
\alias{references-ergm}
\alias{ergm.references}
\alias{references.ergm}
\alias{Poisson}
\alias{Binomial}
\alias{Geometric}
\alias{InitErgmReference.Poisson}
\alias{InitErgmReference.Binomial}
\alias{InitErgmReference.Geometric}
\alias{InitWtErgmProposal.Poisson}
\alias{InitWtErgmProposal.ZIPoisson}
\alias{InitWtErgmProposal.PoissonTNT}
\alias{InitWtErgmProposal.PoissonNonObserved}
\alias{InitWtErgmProposal.Binomial}
\alias{InitWtErgmProposal.BinomialNonObserved}
\alias{InitWtErgmProposal.Geometric}
\alias{InitWtErgmProposal.GeometricNonObserved}

\docType{methods}

\title{Reference Measures for Exponential-Family Random Graph Models for
  Counts}

\description{
  This page describes the possible reference measures (baseline distributions)
  for modeling count data found in the \code{\link[=ergm.count-package]{ergm.count}} package.

  Each of these is specified on the RHS of a one-sided formula passed as
  the \code{reference} argument to \code{\link[ergm]{ergm}}. 
  See the \code{\link[ergm]{ergm}} documentation for a complete
  description of how reference measures are specified.

}

\section{Known issues}{
  
  \subsection{Parameter space constraints}{
    
    Poisson- and geometric-reference ERGMs have an unbouded sample
    space. This means that the parameter space may be constrained in
    complex ways that depend on the terms used in the model. At this
    time \code{\link{ergm}} has no way to detect when a parameter
    configuration had strayed outside of the parameter space, but it may
    be noticeable on a runtime trace plot (activated via
    \code{MCMC.runtime.traceplot} control parameter), when the simulated
    values keep climbing upwards. (See Krivitsky (2012) for a further
    discussion.)

    A possible remedy if this appears to occur is to try lowering the
    control parameter \code{MCMLE.steplength}.

  }

  \subsection{\code{\link[ergm:control.ergm]{MCMLE.trustregion}}}{

    Because Monte Carlo MLE's approximation to the likelihood becomes less
    accurate as the estimate moves away from the one used for the
    sample, \code{\link{ergm}} limits how far the optimization can move
    the estimate for every iteration: the log-likelihood may not change
    by more than \code{\link[ergm:control.ergm]{MCMLE.trustregion}} control
    parameter, which defaults to 20. This is an adequate value for
    binary ERGMs, but because each dyad in a valued ERGM contains more
    information, this number may be too small, resulting in
    unnecessarily many iterations needed to find the MLE.

    Automatically setting \code{\link[ergm:control.ergm]{MCMLE.trustregion}} is
    work in progress, but, in the meantime, you may want to set it to a
    high number (e.g., 1000).
    
  }
  
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
      \eqn{1}.

      Three proposal functions are currently implemented, two of them
      designed to improve mixing for sparse networks. They can can be
      selected via the \code{MCMC.prop.weights=} control parameter. The
      sparse proposals work by proposing a jump to 0. Both of them take
      an optional proposal argument \code{p0} (i.e.,
      \code{MCMC.prop.args=list(p0=...)})  specifying the probability of
      such a jump. However, the way in which they implement it are
      different:

      \describe{
	\item{\code{"random"}}{ Select a dyad \eqn{(i,j)}at random, and draw the
          proposal \eqn{y_{i,j}^\star \sim \mathrm{Poisson}_{\ne
              y_{i,j}}(y_{i,j}+0.5)} (a Poisson distribution with mean
          slightly higher than the current value and conditional on not
          proposing the current value).}
	
        \item{\code{"0inflated"}}{As \code{"random"} but, with
          probability \code{p0}, propose a jump to 0 instead of a
          Poisson jump (if not already at 0). If \code{p0} is not given,
          defaults to the "surplus" of 0s in the observed network,
          relative to Poisson. }
	
	\item{\code{"TNT"} (the default)}{As \code{"0inflated"} but
          instead of selecting a dyad at random, select a tie with
          probability \code{p0}, and a random dyad otherwise, as with
          the binary TNT. Currently, \code{p0} defaults to 0.2.}

      }
    }

    \item{\code{Geometric}}{\emph{Geometric-reference ERGM:} Specifies
      each dyad's baseline distribution to be uniform on the natural
      numbers (and \eqn{0}): \eqn{h(y)=1}. In itself, this
      "distribution" is improper, but in the presence of
      \code{\link[ergm:ergm-terms]{sum}}, a geometric
      distribution is induced. Using \code{\link{CMP}} (in addition to
      \code{\link[ergm:ergm-terms]{sum}}) induces a
      Conway-Maxwell-Poisson distribution that is geometric when its
      coefficient is \eqn{0} and Poisson when its coefficient is
      \eqn{-1}.  }

    \item{\code{Binomial(trials)}}{\emph{Binomial-reference ERGM:}
      Specifies each dyad's baseline distribution to be binomial with
      \code{trials} trials and success probability of \eqn{0.5}:
      \eqn{h(y)=\prod_{i,j}{{\code{trials}}\choose{y_{i,j}}}}.  Using
      \code{\link[ergm:ergm-terms]{valued ERGM terms}} that are
      ``generalized'' from their binary counterparts, with form
      \code{"sum"} (see previous link for the list) produces logistic
      regression.  }

    }
}

\references{
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}

Shmueli G, Minka TP, Kadane JB, Borle S, and Boatwright P (2005). A Useful Distribution for Fitting Discrete Data: Revival of the Conway--Maxwell--Poisson Distribution. \emph{Journal of the Royal Statistical Society: Series C}, 54(1): 127-142.
}

\seealso{ergm, network, \%v\%, \%n\%, sna,
         summary.ergm, print.ergm}
\keyword{models}