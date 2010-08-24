\name{lifetable}
\alias{lifetable}
\title{Construct lifetable from mortality rates}
\usage{
lifetable(data, series = names(data$rate)[1], years = data$year,
                 ages = data$age, max.age = min(100, max(ages)),
                 type = c("period", "cohort"))
}

\arguments{
\item{data}{Demogdata object such as obtained from \code{\link{read.demogdata}}, \code{\link{forecast.fdm}} or \code{\link{forecast.lca}}.
}
\item{series}{Name of series to use.  Default is the first series in data\$rate.}
\item{years}{Vector indicating which years to include in table. Default is all years available.}
\item{ages}{Vector indicating which ages to include in table. Default is all ages available.}
\item{max.age}{Age for last row. Ages beyond this are combined.}
\item{type}{Type of lifetable: \dQuote{period} or \dQuote{cohort}.}
}

\description{Creates lifetable replicating the excel spreadsheet life table. a0 is set to be
\eqn{0.07 + (1.7*m0)}{0.07+1.7m_0} as per Keyfitz.   \eqn{qx = mx/(1 + ((1-ax) *
mx))}{q_x = m_x/(1+[(1-a_x)m_x])} as per Chiang (1984).

\code{lifetable} constructs a lifetable for multiple years.
}

\value{Object of class \dQuote{lifetable} containing the following components:
\item{label}{Name of region from which data are taken.}
\item{series}{Name of series}
\item{age}{Ages for lifetable}
\item{year}{Period years or cohort years}
\item{mx}{Death rate at age x.}
\item{qx}{The probability that an individual of exact age x will die before exact age x+1.}
\item{lx}{Number of survivors to exact age x.  The radix is 1.}
\item{dx}{The number of deaths between exact ages x and x+1.}
\item{Lx}{Number of years lived between exact age x and exact age x+1.}
\item{Tx}{Number of years lived after exact age x.}
\item{ex}{Remaining life expectancy at exact age x.}
}

\seealso{\code{\link{life.expectancy}}}

\author{Heather Booth, Leonie Tickle, Rob J Hyndman, John Maindonald and Timothy Miller}

\references{Chiang CL. (1984) \emph{The life table and its applications}.
    Robert E Krieger Publishing Company: Malabar.

    Keyfitz, N, and Caswell, H. (2005) \emph{Applied mathematical demography}, Springer-Verlag: New York.
}



\examples{
france.lt <- lifetable(fr.mort)
plot(france.lt)

france.LC <- lca(fr.mort)
france.fcast <- forecast(france.LC)
france.lt.f <- lifetable(france.fcast)
plot(france.lt.f,years=2010)
}
\keyword{models}
