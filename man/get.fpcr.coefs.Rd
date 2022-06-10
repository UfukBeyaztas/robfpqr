\name{get.fpcr.coeff}
\alias{get.fpcr.coeff}
\title{Get the estimated regression coefficient functions for scalar-on-function linear and/or quantile regression model}
\description{This function is used to obtain the estimated regression coefficient function \eqn{\beta(t)} or \eqn{\beta(t,~\tau)} for scalar-on-function linear and/or quantile regression model based on output object obtained from \code{\link{fpcr}}).}
\usage{
get.fpcr.coeff(object)
}
\arguments{
  \item{object}{The output object of \code{\link{fpcr}}.}
}
\value{A vector containing the estimated coefficient function \eqn{\beta(t)} or \eqn{\beta(t,~\tau)} depending on \code{model.type}.
}
\author{Ufuk Beyaztas, Mujgan Tez and Han Lin Shang}
\examples{
sim.data <- dgp(n.train = 420, n.test = 180, n.gp = 201, data.type = "contaminated",
            out.type = "yx", out.perc = 0.1)
y <- sim.data$y.train
x <- sim.data$x.train
gp <- seq(0, 1, length.out = 201) # grid points of X
model.linear <- fpcr(y=y, x=x, nbf=20, gp=gp, ncp=4,
                     model.type = "linear")
model.quantile <- fpcr(y=y, x=x, tau=0.5, nbf=20, gp=gp, ncp=4,
                     model.type = "quantile")
linear.coeff.fun <- get.fpcr.coeff(object = model.linear)
quantile.coeff.fun <- get.fpcr.coeff(object = model.quantile)
}
