\name{dgmm}
\alias{dgmm}
\docType{methods}
\title{Density of The Gaussian Mixture Distribution}
\description{
    Density function for the Gaussian mixture distribution.
}
\usage{
    dgmm( x, p, implementation = "C" )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the location of \emph{i}-th component,
        sigma\emph{i} is the scale of \emph{i}-th component.
    }
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{A vector.}
\author{Andrius Merkys}
