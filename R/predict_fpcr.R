predict_fpcr <- function(object, xnew)
{
  gp <- object$model.details[[2]]
  nbf <- object$model.details[[3]]
  PCAcoef <- object$model.details[[4]]
  bs_basis <- object$model.details[[6]]
  coeffs <- object$coeffs

  x.scores.test <- getPCA.test(data = xnew, bs_basis = bs_basis,
                               PCAcoef = PCAcoef, gp = gp)

  preds <- cbind(1,x.scores.test) %*% coeffs

  return(preds)
}
