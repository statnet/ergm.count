#  File R/InitWtErgmTerm.R in package ergm.count, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################
InitWtErgmTerm.CMP<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="sumlogfactorial",
       coef.names="CMP",
       inputs=NULL,
       dependence=FALSE,
       minval=0)
}

InitWtErgmTerm.CMB<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("trials", "coupled"),
                      vartypes = c("numeric", "logical"),
                      defaultvalues = list(NULL, TRUE),
                      required = c(TRUE, FALSE) )
  list(name=if(a$coupled) "CMB" else "CMB2",
       coef.names=if(a$coupled) "CMB" else c("CMB.y","CMB.n-y"),
       inputs=a$trials,
       dependence=FALSE,
       minval=0)
}
