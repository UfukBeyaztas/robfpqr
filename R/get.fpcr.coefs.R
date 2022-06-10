get.fpcr.coeff <- function(object)
{
  gp <- object$model.details$gp
  PCAcoef <- object$model.details$PCAcoef
  evalbase <- object$model.details$evalbase
  coeffs <- as.matrix(object$coeffs[-1])

  f.coeff = evalbase %*% (PCAcoef$coefs %*% coeffs)

  return(f.coeff)
}
