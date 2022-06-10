\name{predict_fpcr}
\alias{predict_fpcr}
\title{Prediction for a scalar-on-function linear and/or quantile regression model based on functional principal component analysis}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted scalar-on-function linear and/or quantile regression model in the output of \code{\link{fpcr}}.}
\usage{
predict_fpcr(object, xnew)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An output object obtained from \code{\link{fpcr}}.}
  \item{xnew}{A matricx consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fpcr}}.}
}
\value{An \eqn{n_{test} \times 1}-dimensional matrix of predicted values of the scalar response variable for the given set of new functional predictor \code{xnew}. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{xnew}.
}
\author{Ufuk Beyaztas, Mujgan Tez and Han Lin Shang}
\examples{
set.seed(123)
sim.data <- dgp(n.train = 420, n.test = 180, n.gp = 201, data.type = "contaminated",
                out.type = "yx", out.perc = 0.1)
y.train <- sim.data$y.train
y.test <- sim.data$y.test
x.train <- sim.data$x.train
x.test <- sim.data$x.test
gp <- seq(0, 1, length.out = 201) # grid points of X
model.linear <- fpcr(y=y.train, x=x.train, nbf=20, gp=gp,
                     ncp=4, model.type = "linear")
model.quantile <- fpcr(y=y.train, x=x.train, tau=0.5, nbf=20,
                       gp=gp, ncp=4, model.type = "quantile")
predict.linear <- predict_fpcr(object=model.linear, xnew=x.test)
predict.quantile <- predict_fpcr(object=model.quantile, xnew=x.test)
# Mean squared prediction error
mspe.linear <- mean((predict.linear-y.test)^2)
mspe.quantile <- mean((predict.quantile-y.test)^2)
round(mspe.linear, 4) # 3.0747
round(mspe.quantile, 4) # 3.3102
}
