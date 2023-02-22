\name{llsmm}
\alias{llsmm}
\docType{methods}
\title{Log-likelihood for Student's t Mixture}
\description{
    Calculates log-likelihood for a given data vector using a Student's t mixture distribution.
}
\usage{
    llsmm( x, p )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 4*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n}, ni1, ni2, ..., ni\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        k\emph{i} is the concentration of \emph{i}-th component and
        ni\emph{i} is the degrees of freedom of \emph{i}-th component.
    }
}
\value{log-likelihood}
\author{Andrius Merkys}
