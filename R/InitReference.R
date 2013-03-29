InitReference.Poisson <- function(lhs.nw, ...){
  list(name="Poisson")
}

InitReference.Binomial <- function(lhs.nw, trials, ...){
  
  list(name="Binomial", trials=trials)
}

InitReference.Geometric <- function(lhs.nw, ...){
  
  list(name="Geometric")
}
