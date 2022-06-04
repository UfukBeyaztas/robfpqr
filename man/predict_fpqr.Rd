\name{predict_fpqr}
\alias{predict_fpqr}
\title{Prediction for a scalar-on-function linear quantile regression model based on functional partial quantile regression}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted scalar-on-function linear quantile regression model in the output of \code{\link{fpqr}}.}
\usage{
predict_fpqr(object, xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{fpqr}}.}
  \item{xnew}{A matricx consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fpqr}}.}
}
\value{An \eqn{n_{test} \times 1}-dimensional matrix of predicted values of the scalar response variable for the given set of new functional predictor \code{xnew}. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{xnew}.}
\author{Ufuk Beyaztas, Mujgan Tez and Han Lin Shang}
\examples{
sim.data <- dgp(n.train = 420, n.test = 180, n.gp = 201, data.type = "contaminated",
            out.type = "yx", out.perc = 0.1)
y.train <- sim.data$y.train
y.test <- sim.data$y.test
x.train <- sim.data$x.train
x.test <- sim.data$x.test
gp <- seq(0, 1, length.out = 201) # grid points of X
model.classic <- fpqr(y=y.train, x=x.train, tau=0.5, h=4, nbasis=20,
                      gp=gp, method.type = "classical")
model.robust <- fpqr(y=y.train, x=x.train,tau=0.5, h=4, nbasis=20, gp=gp,
                     method.type = "robust")
predict.classic <- predict_fpqr(object=model.classic, xnew=x.test)
predict.robust <- predict_fpqr(object=model.robust, xnew=x.test)
}
