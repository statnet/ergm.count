%  File man/ergm-references.Rd in package ergm.count, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2008-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-references}
\alias{ergm-references}
\alias{references-ergm}
\alias{ergm.references}
\alias{references.ergm}
\alias{Poisson}
\alias{Binomial}
\alias{Geometric}
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
}


\section{Possible reference measures to represent baseline
  distributions}{

}

\references{
Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \doi{10.1214/12-EJS696}

Shmueli G, Minka TP, Kadane JB, Borle S, and Boatwright P (2005). A Useful Distribution for Fitting Discrete Data: Revival of the Conway--Maxwell--Poisson Distribution. \emph{Journal of the Royal Statistical Society: Series C}, 54(1): 127-142.
}

\seealso{ergm, network, \%v\%, \%n\%, sna,
         summary.ergm, print.ergm}
\keyword{models}
