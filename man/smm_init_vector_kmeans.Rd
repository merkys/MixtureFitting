\name{smm_init_vector_kmeans}
\alias{smm_init_vector_kmeans}
\docType{methods}
\title{Estimate Student's t Mixture parameters using Expectation Maximisation.}
\description{
    Estimate an initialization vector for Student's t mixture fitting via Expectation Maximisation.
    R implementation of k-means in kmeans() is used to find data point assignment to clusters.
    s_fit_primitive() is then used to estimate component parameters for each cluster.
}
\usage{
    smm_init_vector_kmeans( x, m )
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
