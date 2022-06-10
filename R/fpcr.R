fpcr <- function(y, x, tau, nbf, gp, ncp, model.type = c("linear", "quantile"))
{
  model.type <- match.arg(model.type)

  pca.results <- getPCA(data = x, nbasis = nbf, ncomp = ncp, gp = gp)
  x.scores <- pca.results$PCAscore

  if(model.type == "linear")
  {
    model <- lm(y~x.scores)
  }else if(model.type == "quantile")
  {
    model <- rq(y~x.scores, tau)
  }

  coeffs <- as.matrix(model$coefficients)
  fitteds <- cbind(1, x.scores) %*% coeffs
  resids <- y - fitteds

  model.details <- list()
  model.details$ncp <- ncp
  model.details$gp <- gp
  model.details$nbf <- nbf
  model.details$PCAcoef <- pca.results$PCAcoef
  model.details$evalbase <- pca.results$evalbase
  model.details$bs_basis <- pca.results$bs_basis



  return(list(y = y, x = x, fitted.values = fitteds, residuals = resids,
              coeffs = coeffs, model.details = model.details))
}
