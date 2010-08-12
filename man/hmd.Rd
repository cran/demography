\name{hmd.mx}
\alias{hmd}
\alias{hmd.mx}
\title{Read Mx demographic data from the Human Mortality Database}
\usage{hmd.mx(country, username, password, label = country)  
}

\arguments{
\item{country}{Directory abbreviation from the HMD. For instance, Australia = "AUS".}
\item{username}{HMD username}
\item{password}{HMD password}
\item{label}{Character string giving name of country from which the data are taken.}
}

\description{Reads "Mx" (1x1) data from the Human Mortality Database (\url{http://www.mortality.org}) and constructs a demogdata object
suitable for plotting using \code{\link{plot.demogdata}} and fitting an LC or BMS model using \code{\link{lca}} or an 
FDA model using \code{\link{fdm}}.}

\details{In order to read the data, users are required to create their account via the HMD website (\url{http://www.mortality.org}), 
and obtain a valid username and password.}

\value{Object of class \dQuote{demogdata} with the following components:
\item{year}{Vector of years}
\item{age}{Vector of ages}
\item{rate}{A list containing one or more rate matrices with one age group per row and one column per year.}
\item{pop}{A list of the same form as \code{rate} but containing population numbers instead of demographic rates.}
\item{type}{Type of object: \dQuote{mortality}, \dQuote{fertility} or \dQuote{migration}.}
\item{label}{label}
}

\seealso{\code{\link{demogdata}},\code{\link{read.demogdata}}}

\examples{
\dontrun{norway <- hmd.mx("NOR", username, password, "Norway")
summary(norway)
}
}

\author{Rob J Hyndman}
\keyword{manip}
