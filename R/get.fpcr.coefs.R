get.fpcr.coefs <- function(object)
{
  gp <- object$model.details[[2]]
  PCAcoef <- object$model.details[[4]]
  evalbase <- object$model.details[[5]]
  coeffs <- as.matrix(object$coeffs[-1])

  f.coeff = evalbase %*% (PCAcoef$coefs %*% coeffs)

  return(f.coeff)
}
