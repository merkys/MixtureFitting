\name{polyroot_NR}
\alias{polyroot_NR}
\docType{methods}
\title{Find one real polynomial root using Newton--Raphson method.}
\description{
    Finds one real polynomial root using Newton--Raphson method,
    implemented according to Wikipedia.
}
\usage{
    polyroot_NR( p, init = 0, epsilon = 1e-6, debug = FALSE, implementation = "C" )
}
\arguments{
    \item{p}{vector of polynomial coefficients.}
    \item{init}{initial value.}
    \item{epsilon}{tolerance threshold for convergence.}
    \item{debug}{flag to turn the debug prints on/off.}
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{
    Real polynomial root.
}
\references{
    Users of Wikipedia.
    Newton's method.
    \url{https://en.wikipedia.org/w/index.php?title=Newton\%27s_method&oldid=710342140}
}
\author{Andrius Merkys}
