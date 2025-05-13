\name{smm_fit_em}
\alias{smm_fit_em}
\docType{methods}
\title{Estimate Student's t Mixture parameters using Expectation Maximization.}
\description{
    Estimates parameters for Student's t mixture using Expectation Maximization algorithm.
    Calls smm_fit_em_APK10().
}
\usage{
    smm_fit_em( x, p, ... )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialization vector of 4*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n}, ni1, ni2, ..., ni\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        k\emph{i} is the concentration of \emph{i}-th component and
        ni\emph{i} is the degrees of freedom of \emph{i}-th component.
    }
    \item{...}{
        additional arguments passed to smm_fit_em_GNL08().
    }
}
\value{
    Vector of mixture parameters, whose structure is the same as of input
    parameter's p.
}
\author{Andrius Merkys}
