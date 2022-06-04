\name{fpcr}
\alias{fpcr}
\title{Scalar-on-function linear quantile regression based on principal component analysis}
\description{This function can be used to perform both scalar-on-function linear regression model \deqn{
Y = \int X(t) \beta(t) dt + \epsilon,} and scalar-on-function linear quantile regression model \deqn{
Q_{\tau}(Y|X) = \int X(t) \beta_{\tau}(t) dt} based on the functional principal component decomposition of the functional predictor.}
\usage{fpcr(y, x, tau, nbf, gp, ncp, model.type = c("linear", "quantile"))}
\arguments{
  \item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix with dimension \eqn{n \times p} containing the observations of functional predictor \eqn{X(t)}, where \eqn{n} is the sample size and \eqn{p} denotes the number of grid points for \eqn{X(t)}.}
  \item{tau}{Quantile level.}
  \item{nbf}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional principal components for the functional predictor \eqn{X(t)}.}
  \item{gp}{A vector containing the grid points of the functional predictor \eqn{X(t)}.}
  \item{ncp}{A numeric value denoting the number of functional principal components to be computed for the functional predictor \eqn{X(t)}.}
  \item{model.type}{Fitting model used to estimate the scalar-on-function regression model. Possibilities are "linear" and "quantile".}
}
\value{A list object with the following components:
\item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}.}
\item{x}{A matrix with dimension \eqn{n \times p} containing the observations of functional predictor \eqn{X(t)}.}
\item{fitted.values}{An \eqn{n \times 1}-dimensional matrix containing the fitted values of the scalar response.}
\item{residuals}{An \eqn{n \times 1}-dimensional matrix containing the residuals.}
\item{coeffs}{A vector containing the estimate of parameters of the regression model coducted between the scalar response and principal component scores of the functional predictor.}
\item{model.details}{A list object containing model details, such as number of basis functions, number of principal components, and grid points used for the functional predictor variable.}
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
}
