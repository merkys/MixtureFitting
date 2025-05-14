\name{gmm_fit_hwhm_spline_deriv}
\alias{gmm_fit_hwhm_spline_deriv}
\docType{methods}
\title{Estimate Gaussian Mixture Parameters Using Derivatives and Half-Width at Half-Maximum Method.}
\description{
    Estimate Gaussian mixture parameters using derivatives and half-width at half-maximum (HWHM) method.
    The method smooths the histogram before attempting to locate the modes.
    Then it describes them using HWHM.
}
\usage{
    gmm_fit_hwhm_spline_deriv( x, y )
}
\arguments{
    \item{x}{data vector}
    \item{y}{response vector for \emph{x}}
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
