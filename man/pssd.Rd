\name{pssd}
\alias{pssd}
\docType{methods}
\title{Penalized Sum of Squared Differences Using Gaussian Mixture Distribution}
\description{
    Given two vectors of same length and a Gaussian mixture, calculate the penalized sum of squared differences (SSD) between the first vector and Gaussian mixture densities measured at points from second vector.
    Penalties are included for proportions and scales that are less than or equal to 0.
}
\usage{
    pssd( x, y, p )
}
\arguments{
    \item{x}{data vector}
    \item{y}{response vector}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the location of \emph{i}-th component,
        sigma\emph{i} is the scale of \emph{i}-th component.
    }
}
\value{Penalized sum of squared differences.}
\author{Andrius Merkys}
