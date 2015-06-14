.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.count",c("statnet"),FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
    packageStartupMessage(paste(c(strwrap(paste("NOTE: The form of the term ",sQuote("CMP")," has been changed in version 3.2 of ",sQuote("ergm.count"),". See the news or help('CMP') for more information.",sep="")),""),collapse="\n"))
  }

  .RegisterMHPs()
  .RegisterConstraintImplications()
  .RegisterInitMethods()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "Poisson", "",  0, "random", "Poisson")
  ergm.MHP.table("c", "Poisson", "",  1, "TNT", "PoissonTNT")
  ergm.MHP.table("c", "Poisson", "",  0, "0inflated", "ZIPoisson")
  ergm.MHP.table("c", "Poisson", "observed",  0, "random", "PoissonNonObserved")

  ergm.MHP.table("c", "Geometric", "",  0, "random", "Geometric")
  ergm.MHP.table("c", "Geometric", "observed",  0, "random", "GeometricNonObserved")

  ergm.MHP.table("c", "Binomial", "",  0, "random", "Binomial")
  ergm.MHP.table("c", "Binomial", "observed",  0, "random", "BinomialNonObserved")
}

.RegisterConstraintImplications <- function(){
}

.RegisterInitMethods <- function(){
  ergm.init.methods("Poisson", c("CD","zeros"))
  ergm.init.methods("Binomial", c("CD","zeros"))
  ergm.init.methods("Geometric", c("CD","zeros"))
}
