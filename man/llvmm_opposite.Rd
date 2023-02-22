\name{llvmm_opposite}
\alias{llvmm_opposite}
\docType{methods}
\title{Opposite Log-likelihood for von Mises Mixture}
\description{
    Calculates opposite log-likelihood for a given data vector using a von Mises mixture distribution.
}
\usage{
    llvmm_opposite( x, p )
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
\value{opposite log-likelihood (negated log-likelihood value)}
\author{Andrius Merkys}
