\name{dgmm}
\alias{dgmm}
\docType{methods}
\title{The Gaussian Mixture Distribution}
\description{
    Density function for the Gaussian mixture distribution.
}
\usage{
    dgmm( x, p, normalise_proportions = FALSE, restrict_sigmas = FALSE,
          implementation = "C" )
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
    \item{normalise_proportions}{
        if TRUE, make component proportions sum up to 1 by dividing each one of them by their sum (R implementation only).
    }
    \item{restrict_sigmas}{
        if TRUE, skip components with scales less or equal to zero (R implementation only).
    }
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{A vector.}
\author{Andrius Merkys}
