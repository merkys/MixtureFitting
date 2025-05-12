\name{gmm_fit_hwhm}
\alias{gmm_fit_hwhm}
\docType{methods}
\title{Estimate Gaussian Mixture Parameters Using Half-Width at Half-Maximum Method.}
\description{
    Estimate Gaussian mixture parameters using half-width at half-maximum (HWHM) method.
    Given a histogram, the method attempts to locate most prominent modes and describe them using HWHM.
}
\usage{
    gmm_fit_hwhm( x, y, n )
}
\arguments{
    \item{x}{data vector}
    \item{y}{response vector for \emph{x}}
    \item{n}{number of mixture components}
}
\value{
    Parameter vector of 3*\emph{n} parameters, where \emph{n} is
    number of mixture components. Structure of p vector is
    p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
    where A\emph{i} is the proportion of \emph{i}-th component,
    mu\emph{i} is the location of \emph{i}-th component,
    sigma\emph{i} is the scale of \emph{i}-th component.
}
\author{Andrius Merkys}
