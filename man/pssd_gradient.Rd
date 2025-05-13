\name{pssd_gradient}
\alias{pssd_gradient}
\docType{methods}
\title{Penalized Sum of Squared Differences Using Gaussian Mixture Distribution}
\description{
    Gradient (derivative) function of pssd().
}
\usage{
    pssd_gradient( x, y, p )
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
\value{Gradient values measured at \emph{x}.}
\author{Andrius Merkys}
