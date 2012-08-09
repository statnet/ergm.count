\name{ergm.count-package}
\alias{ergm.count-package}
\docType{package}
\title{
Fit, Simulate and Diagnose Exponential-Family Models for Networks with
Count Edges 
}
\description{
\code{\link[=ergm.count-package]{ergm.count}} is a set of extensions to
\R package \code{\link[=ergm-package]{ergm}} to fit and simulate from exponential-family random 
graph models for networks whose edge weights are counts. For a list of
functions type \code{help(package='ergm')} \code{help(package='ergm.count')}

For a complete list of the functions, use \code{library(help="ergm")} or \code{library(help="ergm.count")} 
or read the rest of the manual. 

When publishing results obtained using this package, please cite the
original authors as described in \code{citation(package="ergm.count")}.

All programs derived from this package must cite it.

}
\details{
This package contains functions specific to using \code{\link{ergm}} to
model networks whose dyad values are counts. Examples include counts of
conversations, messages, and other interactions.

In particular, this package implements the Poisson reference measure for
use by \code{\link{ergm}} and \code{\link{simulate.ergm}} to fit models
from this family, as well as statistics specific to modeling counts,
such as the \code{\link{CMP}} for the Conway-Maxwell-Poisson Distribution.

For detailed information on how to download and install the software,
go to the \code{\link{ergm}} website:
\url{statnet.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
Pavel N. Krivitsky \email{krivitsky@stat.psu.edu}
}
\references{
Handcock MS, Hunter DR, Butts CT, Goodreau SG, Krivitsky PN and Morris M
(2012). _Fit, Simulate and Diagnose Exponential-Family Models for
Networks_. Version 3.1.  Project home
page at <URL: http://www.statnet.org>, <URL:
CRAN.R-project.org/package=ergm>.

Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}

}
\keyword{ package }
\keyword{ models }
