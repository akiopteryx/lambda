sim.coord <- function(n.obj, n.landmark, var=1, cov=0.5, dim) {
  library(geomorph)
  library(MASS)
  vcm <- matrix(cov, nrow=n.landmark*dim, ncol=n.landmark*dim)
  diag(vcm) <- var
  sim.data <- mvrnorm(n.obj, mu=rep(0, n.landmark*dim), Sigma=vcm)
  sim.data <- arrayspecs(sim.data, n.landmark, dim)
  return(sim.data)
}