% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitWtErgmProposal.R
\name{DiscTNT-ergmProposal}
\alias{DiscTNT-ergmProposal}
\alias{InitWtErgmProposal.DiscTNT}
\title{TNT sampling for some discrete-reference ERGMs}
\description{
This proposal implements \link{Poisson-ergmReference},
\link{Geometric-ergmReference}, \link{Binomial-ergmReference}, and
\link{DiscUnif-ergmReference} when the range of values includes 0,
falling back to \link{Disc-ergmProposal} otherwise, all with arbitrary
dyad-level constraints.
}
\details{
\if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsHtml(ergm:::.buildProposalsList(proposal="DiscTNT"))}}
\if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsText(ergm:::.buildProposalsList(proposal="DiscTNT"))}}
\if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsLatex(ergm:::.buildProposalsList(proposal="DiscTNT"))}}
}
\seealso{
\code{\link{ergmProposal}} for index of proposals currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmProposal", "DiscTNT", "subsection")}

\link{TNT-ergmProposal}
}
\concept{bipartite}
\concept{directed}
\concept{discrete}
\concept{undirected}
\concept{valued}
\keyword{internal}
