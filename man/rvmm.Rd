\name{rvmm}
\alias{rvmm}
\docType{methods}
\title{Random Sample of the von Mises Mixture Model.}
\description{
    Generates a random sample of the von Mises Mixture Model.
}
\usage{
    rvmm( n, p )
}
\arguments{
    \item{n}{sample size}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        k\emph{i} is the concentration of \emph{i}-th component.
    }
}
\value{A vector.}
\references{
    Best & Fisher.
    Efficient Simulation of the von Mises Distribution.
    Journal of the RSS, Series C, 1979, 28, 152-157.
}
\author{Andrius Merkys}
