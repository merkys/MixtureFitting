\name{dvmm}
\alias{dvmm}
\docType{methods}
\title{Density of The von Mises Mixture Model.}
\description{
    Density function for the von Mises Mixture Model.
}
\usage{
    dvmm( x, p )
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
\value{A vector.}
\author{Andrius Merkys}
