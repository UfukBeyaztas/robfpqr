predict_fpcr <- function(object, xnew)
{
  gp <- object$model.details$gp
  nbf <- object$model.details$nbf
  PCAcoef <- object$model.details$PCAcoef
  bs_basis <- object$model.details$bs_basis
  coeffs <- object$coeffs

  x.scores.test <- getPCA.test(data = xnew, bs_basis = bs_basis,
                               PCAcoef = PCAcoef, gp = gp)

  preds <- cbind(1,x.scores.test) %*% coeffs

  return(preds)
}
