dgp <- function(n.train, n.test, n.gp, data.type = c("normal", "contaminated"),
                out.type = c("y","yx"), out.perc)
{
  gpX <- seq(0, 1, length.out = n.gp)
  if(data.type == "contaminated"){
    n.out <- round(n.train * out.perc)
    out.indx <- sample(1:n.train, n.out)
  }

  ksi <- list()
  for(i in 1:5){
    ksi[[i]] <- rnorm((n.train+n.test), 0, sd = (4*i^(-3/2)))
  }

  phi <- list()
  phi2 <- list()
  for(i in 1:5){
    phi[[i]] <- sin(i * pi * gpX) - cos(i * pi * gpX)
    phi2[[i]] <- 2*sin(2*i * pi * gpX) - cos(i * pi * gpX)
  }

  fX <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
  fX2 <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi2[[k]])}))
  if(data.type == "contaminated"){
    if(out.type == "yx"){
      fX[out.indx,] = fX[out.indx,] +
        r_ou(n=n.out, t = gpX, mu=15, alpha = 1, sigma = 1,
             x0=rnorm(n=n.out, mean = 15, sd = 1/sqrt(2*1)))$data
    }
  }
  fXd <- fX

  vBeta0 <- vBeta <- 2* sin(2*pi* gpX)

  fX = fdata(fX, argvals = gpX)
  vBeta = fdata(vBeta, argvals = gpX)

  err = rnorm((n.train+n.test), mean=0, sd=1)
  if(data.type == "contaminated"){
    err[out.indx] <- rnorm(n.out, mean=10, sd=1)
  }

  fY = inprod.fdata(fX, vBeta)
  Ye = fY + err

  x.train <- fXd[1:n.train,]
  x.test <- fXd[-(1:n.train),]
  y.train <- as.matrix(Ye[1:n.train])
  y.test <- as.matrix(Ye[-(1:n.train)])

  return(list("y.train" = y.train, "y.test" = y.test, "x.train" = x.train,
              "x.test" = x.test, "f.coef" = vBeta0))
}
