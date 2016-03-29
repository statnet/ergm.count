%  File inst/NEWS.Rd in package ergm.count, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2008-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm.count-package]{ergm.count}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.

\section{Changes in version 3.2.1}{
  
  %% \subsection{NEW FEATURES}{
  %%   \itemize{
  %%   } 
  %% }
  
  %% \subsection{OTHER USER-VISIBLE CHANGES}{
  %%   \itemize{
  %%   }
  %% }

  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %     \item It is not clear that we need a subsection called "C-LEVEL FACILITIES".
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }

  % \subsection{DEPRECATED & DEFUNCT}{
  %   \itemize{
  %   }
  % }

  \subsection{BUG FIXES}{
    \itemize{
      \item Package now works when loaded but not attached.
      \item Test that was failing with \code{\link[ergm-package]{ergm}}
      3.6.0 has been fixed.
    }
  }
}


\section{Changes in version 3.2}{
  
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{ergm.count} now makes use of contrastive divergence in
      \code{\link[=ergm:ergm-package]{ergm}} to obtain starting values for
      the MCMLE.
      
      \item Poisson-reference ERGMs now have TNT-style proposal,
      optimized for sparse networks.
    } 
  }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item The \code{\link{CMP}} term now adds, rather than subtracts
      \eqn{\sum_{(i,j)}\log(y_{i,j}!)}, whereas 3.1.1 and prior versions
      did the opposite. This means that higher
      \eqn{\theta_{\text{CMP}}} represents overdispersion and lower
      \eqn{\theta_{\text{CMP}}} represents underdispersion.
    }
  }
}


\section{Changes in version 3.1.1}{

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item E-mail address and affiliation of the maintainer have been updated.
    }
  }

}


\section{Changes in version 3.1}{

  \subsection{THIS IS AN INITIAL RELEASE}{
    \itemize{
      \item Although this is an initial release, it is depends very heavily on
      the \code{\link[=ergm-package]{ergm}} package version 3.1, and it is
      thus numbered 3.1.
    }
  }
  
  \subsection{NEW FEATURES}{
    \itemize{
      \item A number of reference measures for count data have been
      implemented. See \code{? 'ergm.count-package'} and \code{? 'ergm-references'} for help.
    } 
  }
}
