get.fpqr.coeff <- function(object)
{
  sinp_mat <- object$model.details$sinp_mat
  evalbase <- object$model.details$evalbase
  gp <- object$model.details$gp
  V <- object$V
  coeffs <- as.matrix(object$pqr.coefs[-1])

  f.coeff <- evalbase %*% (t(solve(sinp_mat)) %*% V %*% coeffs)

  return(f.coeff)
}
