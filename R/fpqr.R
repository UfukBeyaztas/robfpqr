fpqr <- function(y, x, tau, h, nbasis, gp, method.type = c("classical","robust"),
                 probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
                 maxit = 1000, conv = 0.01)
{
  method.type <- match.arg(method.type)

  x0 <- x
  y0 <- y

  BS.sol <- getAmat(data = x, nbf = nbasis, gp = gp)
  x <- BS.sol$Amat
  n <- dim(x)[1]
  p <- dim(y)[2]

  xs <- scale_fun(data = x, method.type = method.type)
  ys <- as.matrix(y)

  if(method.type == "classical")
  {
    m.pqr <- pqr(y = ys, x = xs, tau = tau, h = h, method.type = method.type)
  }else if(method.type == "robust")
  {

    wx <- sqrt(apply(xs^2, 1, sum))
    wx <- wx / median(wx)
    zerows <- vector(length=0)

    probct <- qnorm(probp1)
    hampelb <- qnorm(hampelp2)
    hampelr <- qnorm(hampelp3)
    wx[which(wx <= probct)] <- 1
    wx[which(wx > probct & wx <= hampelb)] <- probct / abs(wx[which(wx > probct & wx <= hampelb)])
    wx[which(wx > hampelb & wx <= hampelr)] <- probct * (hampelr - abs(wx[which(wx > hampelb & wx <= hampelr)]))/
      (hampelr - hampelb) * 1 / abs(wx[which(wx > hampelb & wx <= hampelr)])
    wx[which(wx > hampelr)] <- 0

    if(any(wx < 1e-6))
    {
      w0 <- which(wx < 1e-6)
      wx <- replace(wx, list = w0, values = 1e-6)
    }else {
      wx <- wx
    }

    xs.w <- as.data.frame(xs) * sqrt(wx)
    loops <- 1
    rold <- 10^-5
    difference <- 1

    while((difference > conv) && loops < maxit)
    {
      m.pqr <- pqr(y = ys, x = xs.w, tau = tau, h = h, method.type = method.type)
      yp <- m.pqr$fitted.values
      b <- m.pqr$d.coef
      Tpls <- as.data.frame(m.pqr$T) / sqrt(wx)

      dt <- scale_fun(data = Tpls, method.type = method.type)
      wtn <- sqrt(apply(dt^2, 1, sum))
      wtn <- wtn/median(wtn)

      probct <- qnorm(probp1)
      hampelb <- qnorm(hampelp2)
      hampelr <- qnorm(hampelp3)

      probct <- qchisq(probp1, h)
      hampelb <- qchisq(hampelp2, h)
      hampelr <- qchisq(hampelp3, h)
      wte <- wtn
      wte[which(wtn <= probct)] <- 1
      wte[which(wtn > probct & wtn <= hampelb)] <- probct / abs(wtn[which(wtn > probct & wtn <= hampelb)])
      wte[which(wtn > hampelb & wtn <= hampelr)] <- probct * (hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)])) /
        (hampelr - hampelb) * 1 / abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
      wte[which(wtn > hampelr)] <- 0

      difference <- abs(sum(b^2) - rold)/rold
      rold <- sum(b^2)
      wx <- wx * wte

      if(any(wx < 1e-6)){
        w0 <- which(wx < 1e-6)
        wx <- replace(wx, list = w0, values = 1e-6)
        zerows <- unique(c(zerows, as.numeric(names(w0))))
      }

      xs.w <- as.data.frame(xs) * sqrt(wx)

      loops <- loops + 1
    }

    if (difference > conv){
      warning(paste("Convergence was not ashieved. The scaled difference between norms of the coefficient vectors is ",
                    round(difference, digits=4)))
    }

    w <- wx
    w[zerows] <- 0
    wt <- wte
    wt[zerows] <- 0
  }

  P <- m.pqr$P
  R <- m.pqr$R
  W <- m.pqr$W
  T <- m.pqr$T

  V = t(solve(t(T)%*%T)%*%t(T)%*%BS.sol$Amat)

  final.model <- rq(y0~T, tau)
  pqr.coefs <- as.matrix(final.model$coefficients)
  intercept <- pqr.coefs[1]
  final.coefs <- R %*% as.matrix(pqr.coefs[-1])

  fitteds <- cbind(1, T) %*% pqr.coefs
  resids <- y0 - fitteds

  model.details <- list()
  model.details[[1]] <- BS.sol$Amat
  model.details[[2]] <- BS.sol$bs_basis
  model.details[[3]] <- BS.sol$inp_mat
  model.details[[4]] <- BS.sol$sinp_mat
  model.details[[5]] <- BS.sol$evalbase
  model.details[[6]] <- nbasis
  model.details[[7]] <- gp
  model.details[[8]] <- method.type


  return(list(y = y0, x = x0, fitted.values = fitteds, coeffs = final.coefs, intercept = intercept,
              pqr.coefs = pqr.coefs, V = V, model.details = model.details))
}
