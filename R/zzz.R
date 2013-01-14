.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.count",c("statnet"),FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)

  .RegisterMHPs()
  .RegisterConstraintImplications()
  .RegisterInitMethods()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "Poisson", "",  0, "random", "Poisson")
  ergm.MHP.table("c", "Poisson", "",  0, "0inflated", "ZIPoisson")
  ergm.MHP.table("c", "Poisson", "observed",  0, "random", "PoissonNonObserved")

  ergm.MHP.table("c", "Geometric", "",  0, "random", "Geometric")
  ergm.MHP.table("c", "Geometric", "observed",  0, "random", "GeometricNonObserved")

  ergm.MHP.table("c", "Binomial", "",  0, "random", "Binomial")
  ergm.MHP.table("c", "Binomial", "observed",  0, "random", "BinomialNonObserved")

  ergm.MHP.table("c", "DiscUnif", "",  0, "random", "DiscUnif")
  ergm.MHP.table("c", "DiscUnif", "observed",  0, "random", "DiscUnifNonObserved") 
}

.RegisterConstraintImplications <- function(){
}

.RegisterInitMethods <- function(){
  ergm.init.methods("Poisson", c("zeros"))
  ergm.init.methods("Binomial", c("zeros"))
  ergm.init.methods("Geometric", c("zeros"))
  ergm.init.methods("DiscUnif", c("zeros"))
}
