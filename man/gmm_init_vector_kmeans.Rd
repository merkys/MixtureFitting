\name{gmm_init_vector_kmeans}
\alias{gmm_init_vector_kmeans}
\docType{methods}
\title{Estimate Gaussian Mixture parameters using Expectation Maximisation.}
\description{
    Estimate an initialisation vector for Gaussian mixture fitting using k-means.
    R implementation of k-means in kmeans() is used to find data point assignment to clusters.
}
\usage{
    gmm_init_vector_kmeans( x, m )
}
\arguments{
    \item{x}{data vector}
    \item{m}{number of mixture components}
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
