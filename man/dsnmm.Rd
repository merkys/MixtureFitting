\name{dsnmm}
\alias{dsnmm}
\docType{methods}
\title{Density of The Skew-Normal Mixture Model}
\description{
    Density function for the skew-normal mixture model.
}
\usage{
    dsnmm( x, p, implementation = "C" )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        parameter vector of 4*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( omega1, omega2, ..., omega\emph{n}, dzeta1, dzeta2, ..., dzeta\emph{n}, sigma1, sigma2, ..., sigma\emph{n}, lambda1, lambda2, ..., lambda\emph{n} ),
        where omega\emph{i} is the proportion of \emph{i}-th component,
        dzeta\emph{i} is the location of \emph{i}-th component,
        sigma\emph{i} is the scale of \emph{i}-th component and
        lambda\emph{i} is the shape of \emph{i}-th component.
    }
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{A vector.}
\author{Andrius Merkys}
