\name{llcmm}
\alias{llcmm}
\docType{methods}
\title{Log-likelihood for Cauchy Mixture}
\description{
    Calculates log-likelihood for a given data vector using a Cauchy mixture distribution.
}
\usage{
    llcmm( x, p, implementation = "C" )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, gamma1, gamma2, ..., gamma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        gamma\emph{i} is the Cauchy scale of \emph{i}-th component.
    }
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{log-likelihood}
\author{Andrius Merkys}
