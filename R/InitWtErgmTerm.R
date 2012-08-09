InitWtErgmTerm.CMP<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="nsumlogfactorial",
       coef.names="CMP",
       inputs=NULL,
       dependence=FALSE,
       maxval=0)
}
