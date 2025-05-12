\name{cmm_fit_hwhm_spline_derivatives}
\alias{cmm_fit_hwhm_spline_derivatives}
\docType{methods}
\title{Estimate Cauchy Mixture Parameters Using Derivatives and Half-Width at Half-Maximum Method.}
\description{
    Estimate Cauchy mixture parameters using derivatives and half-width at half-maximum (HWHM) method.
    The method smooths the histogram before attempting to locate the modes.
    Then it describes them using HWHM.
}
\usage{
    cmm_fit_hwhm_spline_derivatives( x, y )
}
\arguments{
    \item{x}{data vector}
    \item{y}{response vector for \emph{x}}
}
\value{
    Parameter vector of 3*\emph{n} parameters, where \emph{n} is
    number of mixture components. Structure of p vector is
    p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, gamma1, gamma2, ..., gamma\emph{n} ),
    where A\emph{i} is the proportion of \emph{i}-th component,
    mu\emph{i} is the center of \emph{i}-th component and
    gamma\emph{i} is the Cauchy scale of \emph{i}-th component.
}
\author{Andrius Merkys}
