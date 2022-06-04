\name{dgp}
\alias{dgp}
\title{Generate a dataset for the scalar-on-function regression model}
\description{This function can be used to generate a dataset for the scalar-on-function regression model \deqn{
Y = \int X(t) \beta(t) dt + \epsilon,} where \eqn{Y} denotes the scalar response, \eqn{X(t)} denotes the functional predictor, \eqn{\beta(t)} denotes the regression coefficient function, and \eqn{\epsilon} is the error process.}
\usage{
dgp(n.train, n.test, n.gp, data.type = c("normal", "contaminated"),
                out.type = c("y","yx"), out.perc)
}
\arguments{
  \item{n.train}{An integer, specifying the number of observations for each variable to be generated in the training sample.}
  \item{n.test}{An integer, specifying the number of observations for each variable to be generated in the test sample.}
  \item{n.gp}{An integer, denoting the number of grid points, i.e., a fine grid on the interval [0, 1].}
  \item{data.type}{Data type to be generated. Possibilities are "normal" and "contaminated".}
  \item{out.type}{Outlier type to be generated. Possibilities are "y" and "yx".}
  \item{out.perc}{A numeric value between 0 and 1 specifying the outlier percentage.}
}
\details{In the data generation process, first, the  functional predictor is generated based on the following process: \deqn{
X(t) = \sum_{j=1}^5 \kappa_j v_j(t),} where \eqn{ \kappa_j } is a vector generated from a Normal distribution with mean zero and variance \eqn{2 j^{-3/2}} and \deqn{v_j(t) = \sin(j \pi t) - \cos(j \pi t).} The  regression coefficient function is generated from \eqn{2 \sin(2 \pi t)}. The error process is generated from the standard normal distribution.

If \code{data.type = "normal"}, then the data are generated as above. On the other hand, if \code{data.type = "contaminated"}, then \eqn{n.train \times out.perc} part of the data in the training sample is contaminated by outliers. If \code{out.type = "y"}, then only the scalar response variable is contaminated by outiers. While doing so, \eqn{n.train \times out.perc} of \eqn{\epsilon} is generated from a normal distribution with mean 10 and variance 1. If \code{out.type = "yx"}, then both the scalar response and functional predictor are contaminated by outliers. In doing so, \eqn{n.train \times out.perc}  of \eqn{X(t)} are contaminated an Ornstein-Uhlenbeck process. In addition, \eqn{n.train \times out.perc} of \eqn{\epsilon} is generated from a normal distribution with mean 10 and variance 1. All the functional predictors are generated equally spaced point in the interval \eqn{[0, 1]}.
}
\value{A list object with the following components:
\item{y.train}{An \eqn{n.train \times 1}-dimensional matrix containing the observations of simulated scalar response variable in the training sample.}
\item{y.test}{An \eqn{n.test \times 1}-dimensional matrix containing the observations of simulated scalar response variable in the test sample.}
\item{x.train}{A matrix with dimension \eqn{n.train \times n.gp} containing the observations of simulated functional predictor variable in the training sample.}
\item{x.test}{A matrix with dimension \eqn{n.test \times n.gp} containing the observations of simulated functional predictor variable in the test sample.}
\item{f.coef}{A vector with length \eqn{n.gp} containing the generated regression coefficient function.}
}
\author{Ufuk Beyaztas, Mujgan Tez and Han Li Shang}
\examples{
library(fda.usc)
sim.data <- dgp(n.train = 420, n.test = 180, n.gp = 201, data.type = "contaminated",
            out.type = "yx", out.perc = 0.1)
y.train <- sim.data$y.train
y.test <- sim.data$y.test
x.train <- sim.data$x.train
x.test <- sim.data$x.test
coeffs <- sim.data$f.coef
plot(y.train, type = "p", pch = 16, xlab = "Index", ylab = "", main = "Response")
fx <- fdata(x.train, argvals = seq(0, 1, length.out = 201))
plot(fx, lty = 1, ylab = "", xlab = "Grid point",
     main = expression(X(t)), mgp = c(2, 0.5, 0))
}
