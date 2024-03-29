% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lca.R
\name{lca}
\alias{lca}
\alias{bms}
\title{Model mortality or fertility data using Lee-Carter approach}
\usage{
lca(
  data,
  series = names(data$rate)[1],
  years = data$year,
  ages = data$age,
  max.age = 100,
  adjust = c("dt", "dxt", "e0", "none"),
  chooseperiod = FALSE,
  minperiod = 20,
  breakmethod = c("bai", "bms"),
  scale = FALSE,
  restype = c("logrates", "rates", "deaths"),
  interpolate = FALSE
)

bms(
  data,
  series = names(data$rate)[1],
  years = data$year,
  ages = data$age,
  max.age = 100,
  minperiod = 20,
  breakmethod = c("bms", "bai"),
  scale = FALSE,
  restype = c("logrates", "rates", "deaths"),
  interpolate = FALSE
)
}
\arguments{
\item{data}{demogdata object of type \dQuote{mortality} or
\dQuote{fertility}. Output from read.demogdata.}

\item{series}{name of series within data containing mortality or fertility
values (1x1)}

\item{years}{years to include in fit. Default: all available years.}

\item{ages}{ages to include in fit. Default: all available ages up to
\code{max.age}.}

\item{max.age}{upper age to include in fit. Ages beyond this are collapsed
into the upper age group.}

\item{adjust}{method to use for adjustment of coefficients \eqn{k_t kt}.
Possibilities are \dQuote{dxt} (BMS method), \dQuote{dt} (Lee-Carter
method), \dQuote{e0} (method based on life expectancy) and \dQuote{none}.
Defaults are \dQuote{dxt} for \code{bms()} and \dQuote{dt} for
\code{lca()}.}

\item{chooseperiod}{If TRUE, it will choose the best fitting period.}

\item{minperiod}{Minimum number of years to include in fitting period if
chooseperiod=TRUE.}

\item{breakmethod}{method to use for identifying breakpoints if
chooseperiod=TRUE. Possibilities are \dQuote{bai} (Bai's method computed
using \code{\link[strucchange]{breakpoints}} in the strucchange package)
and \dQuote{bms} (method based on mean deviance ratios described in BMS).}

\item{scale}{If TRUE, it will rescale bx and kt so that kt has drift
parameter = 1.}

\item{restype}{method to use for calculating residuals. Possibilities are
\dQuote{logrates}, \dQuote{rates} and \dQuote{deaths}.}

\item{interpolate}{If TRUE, it will estimate any zero mortality or fertility
rates using the same age group from nearby years.}
}
\value{
Object of class \dQuote{lca} with the following components:
\item{label}{Name of region}
\item{age}{Ages from \code{data} object.}
\item{year}{Years from \code{data} object.}
\item{<series>}{Matrix of mortality or fertility data as contained in \code{data}. It takes the name given by the series argument.}
\item{ax}{Average deathrates across fitting period}
\item{bx}{First principal component in Lee-Carter model}
\item{kt}{Coefficient of first principal component}
\item{residuals}{Functional time series of residuals.}
\item{fitted}{Functional time series containing estimated mortality or fertility rates from model}
\item{varprop}{Proportion of variance explained by model.}
\item{y}{The data stored as a functional time series object.}
\item{mdev}{Mean deviance of total and base lack of fit, as described in Booth, Maindonald and Smith.}
}
\description{
Lee-Carter model of mortality or fertility rates. \code{lca} produces a
standard Lee-Carter model by default, although many other options are
available. \code{bms} is a wrapper for \code{lca} and returns a model based
on the Booth-Maindonald-Smith methodology.
}
\details{
All mortality or fertility data are assumed to be in matrices of
mortality or fertility rates within \code{data$rate}. Each row is one age group
(assumed to be single years). Each column is one year. The
function produces a model for the \code{series} mortality or fertility rate matrix
within \code{data$rate}. Forecasts from this model can be obtained using \code{\link{forecast.lca}}.
}
\examples{
\dontrun{
france.LC1 <- lca(fr.mort, adjust="e0")
plot(france.LC1)
par(mfrow=c(1,2))
plot(fr.mort,years=1953:2002,ylim=c(-11,1))
plot(forecast(france.LC1,jumpchoice="actual"),ylim=c(-11,1))

france.bms <- bms(fr.mort, breakmethod="bai")
fcast.bms <- forecast(france.bms)
par(mfrow=c(1,1))
plot(fcast.bms$kt)
}
}
\references{
Booth, H., Maindonald, J., and Smith, L. (2002) Applying Lee-Carter
under conditions of variable mortality decline. \emph{Population Studies}, \bold{56}, 325-336.

Lee, R.D., and Carter, L.R. (1992) Modeling and forecasting US mortality. \emph{Journal of
  the American Statistical Association}, \bold{87}, 659-671.
}
\seealso{
\code{\link{forecast.lca}}, \code{\link{plot.lca}}, \code{\link{summary.lca}}, \code{\link{fdm}}
}
\author{
Heather Booth, Leonie Tickle, John Maindonald and Rob J Hyndman.
}
\keyword{models}
