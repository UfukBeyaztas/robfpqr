predict_fpqr <- function(object, xnew)
{
  nbasis <- object$model.details[[6]]
  gp <- object$model.details[[7]]
  method.type <- object$model.details[[8]]
  coeffs <- object$coeffs
  intercept <- object$intercept
  
  BS.sol.test <- getAmat(data = xnew, nbf = nbasis, gp = gp)
  Amat.test <- BS.sol.test$Amat
  Amat.test <- scale_fun(data = Amat.test, method.type = method.type)
  predictions <- Amat.test %*% coeffs + intercept
  
  return(predictions)
}