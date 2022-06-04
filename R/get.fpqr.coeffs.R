get.fpqr.coeffs <- function(object)
{
  sinp_mat <- object$model.details[[4]]
  evalbase <- object$model.details[[5]]
  gp <- object$model.details[[7]]
  V <- object$V
  coeffs <- as.matrix(object$pqr.coefs[-1])
  
  f.coeff <- evalbase %*% (solve(sinp_mat) %*% V %*% coeffs)
  
  return(f.coeff)
}