\name{llgmm_opposite}
\alias{llgmm_opposite}
\docType{methods}
\title{Opposite Log-likelihood for Gaussian Mixture}
\description{
    Calculates opposite log-likelihood for a given data vector using a Gaussian mixture distribution.
}
\usage{
    llgmm_opposite( x, p )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        sigma\emph{i} is the scale of \emph{i}-th component.
    }
}
\value{opposite log-likelihood (negated log-likelihood value)}
\author{Andrius Merkys}
