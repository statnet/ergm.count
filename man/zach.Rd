%  File man/zach.Rd in package ergm.count, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{zach}
\alias{zach}
\docType{data}
\title{
  Karate club social network of Zachary (1977)
}
\description{

  Zachary (1977) reported observations of social relations in a
university karate club, with membership that varied between 50 and 100,
of whom 34 individuals: 32 ordinary club members and officers, the club
president ("John A."), and the part-time instructor ("Mr. Hi");
consistently interacted outside of the club. Over the course of the
study, the club divided into two factions, and, ultimately, split into
two clubs, one led by Hi and the other by John and the original club's
officers. The split was driven by a disagreement over whether Hi could
unilaterally change the level of compensation for his services.

Zachary identifies the faction with which each of the 34 actors was aligned
and how strongly and reports, for each pair of actors, the count of social
contexts in which they interacted. The 8 contexts recorded were
\itemize{
  \item academic classes at the university;
  \item Hi's private karate studio in his night classes;
  \item Hi's private karate studio where he taught on weekends;
  \item student-teaching at Hi's studio;
  \item the university rathskeller (bar) located near the karate club;
  \item a bar located near the university campus;
  \item open karate tournaments in the area; and
  \item intercollegiate karate tournaments.
}

The highest number of contexts of interaction for a pair of individuals
that was observed was 7.
}
\usage{data(zach)}
\format{
  The data are represented as a \code{\link{network}} object, with an edge
attribute \code{contexts}, giving the number of contexts of interaction
for that pair of actors. In addition, the following vertex attributes
are provided:
\describe{
  \item{\code{club}:}{the club in which the actor ended up;}
  \item{\code{faction}:}{faction alignment of the actor as recorded by Zachary} 
  \item{\code{faction.id}}{faction alignment coded numerically, as 
    \eqn{-2} (strongly Mr. Hi's), \eqn{-1} (weakly Mr. Hi's), \eqn{0}
    (neutral), \eqn{+1} (weakly John's), and \eqn{+2} (strongly John's);} 
  \item{\code{role}}{role of the actor in the network (Instructor,
    Member, or President)}
}
}
\source{
  Zachary, WW (1977). An Information Flow Model for Conflict and Fission
  in Small Groups. \emph{Journal of Anthropological Research}, 33(4),
  452-473.

  Sociomatrix in machine-readable format was retrieved from \url{http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm}.
}

\references{
  Zachary, WW (1977). An Information Flow Model for Conflict and Fission
  in Small Groups. \emph{Journal of Anthropological Research}, 33(4),
  452-473.
}

\examples{
\donttest{
data(zach)

# Fit a binomial-reference ERGM.

zach.fit1 <- ergm(zach~nonzero+sum+nodefactor("role",base=2)+absdiffcat("faction.id"),
                  response="contexts", reference=~Binomial(8),
                  control=control.ergm(MCMLE.trustregion=1000))

mcmc.diagnostics(zach.fit1)

summary(zach.fit1)

\dontrun{
# This is much slower.
zach.fit2 <- ergm(zach~nonzero+sum+nodefactor("role",base=2)+transitiveties,
                  response="contexts", reference=~Binomial(8),
                  control=control.ergm(MCMLE.trustregion=1000))

mcmc.diagnostics(zach.fit2)

summary(zach.fit2)
}
}
}
\keyword{datasets}
