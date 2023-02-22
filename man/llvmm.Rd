\name{llvmm}
\alias{llvmm}
\docType{methods}
\title{Log-likelihood for von Mises Mixture}
\description{
    Calculates log-likelihood for a given data vector using a von Mises mixture distribution.
}
\usage{
    llvmm( x, p )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        k\emph{i} is the concentration of \emph{i}-th component.
    }
}
\value{log-likelihood}
\author{Andrius Merkys}
