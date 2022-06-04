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
  model.details[[1]] <- ncp
  model.details[[2]] <- gp
  model.details[[3]] <- nbf
  model.details[[4]] <- pca.results$PCAcoef
  model.details[[5]] <- pca.results$evalbase
  model.details[[6]] <- pca.results$bs_basis

  return(list(y = y, x = x, fitted.values = fitteds, residuals = resids,
              coeffs = coeffs, model.details = model.details))
}
