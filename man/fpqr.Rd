\name{fpqr}
\alias{fpqr}
\title{Functional partial quantile regression}
\description{This function is used to perform scalar-on-function linear quantile regression model \deqn{
Q_{\tau}(Y|X) = \int X(t) \beta_{\tau}(t) dt} based on the partial quantile regression.}
\usage{
fpqr(y, x, tau, h, nbasis, gp, method.type = c("classical","robust"),
     probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
     maxit = 1000, conv = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix with dimension \eqn{n \times p} containing the observations of functional predictor \eqn{X(t)}, where \eqn{n} is the sample size and \eqn{p} denotes the number of grid points for \eqn{X(t)}.}
  \item{tau}{Quantile level.}
  \item{h}{A numeric value denoting the number of functional partial quantile regression components to be computed.}
  \item{nbasis}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional partial quantile regression components.}
  \item{gp}{A vector containing the grid points of the functional predictor \eqn{X(t)}.}
  \item{method.type}{Method type used to estimate the scalar-on-function linear quantile regression model. Possibilities are "classical" and "robust".}
  \item{probp1}{A numeric value used to determine the first outlier cutoff point for the weights.}
  \item{hampelp2}{A numeric value used to determine the first outlier cutoff point for the weights.}
  \item{hampelp3}{A numeric value used to determine the third outlier cutoff point for the weights.}
  \item{maxit}{An integer value defining the maximum iteration used to achieve convergence.}
  \item{conv}{A numeric value used for the precision of the coefficient estimate.}
}
\details{If \code{method.type = "classical"}, then, the partial quantile regression algorithm of Dodge and Whittaker (2009) is used to obtain functional partial quantile regression components.

If \code{method.type = "robust"}, then, the partial quantile regression algorithm of Dodge and Whittaker (2009) along with a modified version of the iteratively reweighting algorithm of Serneels et al. (2005) is used to obtain functional partial quantile regression components.}
\value{A list object with the following components:
\item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}.}
\item{x}{A matrix with dimension \eqn{n \times p} containing the observations of functional predictor \eqn{X(t)}.}
\item{fitted.values}{An \eqn{n \times 1}-dimensional matrix containing the fitted values of the scalar response.}
\item{residuals}{An \eqn{n \times 1}-dimensional matrix containing the residuals.}
\item{coeffs}{A vector containing the estimate of parameters of the regression model coducted between the scalar response and principal component scores of the functional predictor.}
\item{intercept}{A numeric value containing the estimated intercept parameter.}
\item{pqr.coefs}{A vector containing the estimated regression parameter for the regression problem of scalar response on the partial quantile regression components.}
\item{V}{A matrix whose rows are the eigenvectors}
\item{model.details}{A list object containing model details, such as number of basis functions, number of partial quantile components, and grid points used for the functional predictor variable.}
}
\references{
Y. Dodge and J. Whittaker (2009), "Partial quantile regression, \emph{Metrika}, \bold{70}(1), 35–57.

S. Serneels and C. Croux and P. Filzmoser and P. J. V. Espen (2005), "Partial robust M-regression", \emph{Chemometrics and Intelligent Laboratory Systems}, \bold{79}(1-2), 55–64.
}
\author{Ufuk Beyaztas, Mujgan Tez and Han Lin Shang}
\examples{
sim.data <- dgp(n.train = 420, n.test = 180, n.gp = 201, data.type = "contaminated",
            out.type = "yx", out.perc = 0.1)
y <- sim.data$y.train
x <- sim.data$x.train
gp <- seq(0, 1, length.out = 201) # grid points of X
model.classic <- fpqr(y=y, x=x, tau=0.5, h=4, nbasis=20, gp=gp, method.type = "classical")
model.robust <- fpqr(y=y, x=x,tau=0.5, h=4, nbasis=20, gp=gp, method.type = "robust")
}
