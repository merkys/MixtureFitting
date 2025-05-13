\name{gmm_init_vector_quantile}
\alias{gmm_init_vector_quantile}
\docType{methods}
\title{Estimate Gaussian Mixture parameters using Expectation Maximisation.}
\description{
    Estimate an initialization vector for Gaussian mixture fitting using (weighted) quantiles.
    Proportions and scales are set to equal, centers are placed at equispaced quantiles.
}
\usage{
    gmm_init_vector_quantile( x, m, w = numeric() )
}
\arguments{
    \item{x}{data vector}
    \item{m}{number of mixture components}
    \item{w}{weight vector}
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
