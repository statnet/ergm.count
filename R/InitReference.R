InitReference.Poisson <- function(lhs.nw, ...){
  list(name="Poisson")
}

InitReference.Binomial <- function(lhs.nw, trials, ...){
  
  list(name="Binomial", trials=trials)
}

InitReference.Geometric <- function(lhs.nw, ...){
  
  list(name="Geometric")
}

InitReference.DiscUnif <- function(lhs.nw, a, b, ...){
  
  list(name="DiscUnif", a=a, b=b)
}
