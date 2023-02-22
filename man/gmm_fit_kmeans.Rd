\name{gmm_fit_kmeans}
\alias{gmm_fit_kmeans}
\docType{methods}
\title{Estimate Gaussian Mixture parameters from kmeans.}
\description{
    Estimates parameters for Gaussian mixture using kmeans.
}
\usage{
    gmm_fit_kmeans( x, n )
}
\arguments{
    \item{x}{data vector}
    \item{n}{number of mixture components}
}
\value{
    Vector of 3*\emph{n} mixture parameters, where \emph{n} is
    number of mixture components. Structure of p vector is
    p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
    where A\emph{i} is the proportion of \emph{i}-th component,
    mu\emph{i} is the location of \emph{i}-th component,
    sigma\emph{i} is the scale of \emph{i}-th component.
}
\author{Andrius Merkys}
