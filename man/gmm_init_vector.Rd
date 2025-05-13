\name{gmm_init_vector}
\alias{gmm_init_vector}
\docType{methods}
\title{Estimate Gaussian Mixture parameters using Expectation Maximisation.}
\description{
    Estimate an initialization vector for Gaussian mixture fitting via Expectation Maximisation.
    Proportions and scales are set to equal, centers are equispaced through the whole domain of input sample.
}
\usage{
    gmm_init_vector( x, n, implementation = "C" )
}
\arguments{
    \item{x}{data vector}
    \item{n}{number of mixture components}
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
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
