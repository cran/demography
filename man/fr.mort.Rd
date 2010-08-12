\name{fr.mort}
\alias{fr.mort}

\docType{data}

\title{French mortality data}

\description{Age-specific mortality rates and population for France.}

\details{French mortality rates and populations (1899-2005) for ages 0-110. Data taken from the Human Mortality Database
on 20 February 2008.
}

\format{Object of class \code{demogdata} containing the following components:
\describe{
\item{year}{Vector of years}
\item{age}{Vector of ages}
\item{rate}{List of matrices containing rates with with one age group per row and one column per year.
    Matrices: \code{total}, \code{female}, \code{male}.}
\item{pop}{Population data in same form as \code{rate}.}
\item{type}{Type of object. In this case, \dQuote{mortality}.}
\item{label}{Character string giving area from which data are taken. In this case, \dQuote{France}.}
}}

\source{The Human Mortality Database (\url{http://www.mortality.org}).}

\author{Rob J Hyndman}

\examples{
plot(fr.mort,years=1950:1997)
}
\keyword{data}
