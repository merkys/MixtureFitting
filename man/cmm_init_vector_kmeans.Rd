\name{cmm_init_vector_kmeans}
\alias{cmm_init_vector_kmeans}
\docType{methods}
\title{Estimate Cauchy Mixture parameters using Expectation Maximization.}
\description{
    Estimate an initialization vector for Cauchy mixture fitting using k-means.
    R implementation of k-means in kmeans() is used to find data point assignment to clusters.
    Then several iterations of Cauchy mixture fitting (per Nahy 2006) is used to derive mixture parameters.
}
\usage{
    cmm_init_vector_kmeans( x, m, iter.cauchy = 20 )
}
\arguments{
    \item{x}{data vector}
    \item{m}{number of mixture components}
    \item{iter.cauchy}{number of iterations to fit a single Cauchy component.}
}
\value{
    Parameter vector of 3*\emph{n} parameters, where \emph{n} is
    number of mixture components. Structure of p vector is
    p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, gamma1, gamma2, ..., gamma\emph{n} ),
    where A\emph{i} is the proportion of \emph{i}-th component,
    mu\emph{i} is the center of \emph{i}-th component and
    gamma\emph{i} is the Cauchy scale of \emph{i}-th component.
}
\references{
    Ferenc Nahy. Parameter Estimation of the Cauchy Distribution in
    Information Theory Approach (2006).
    Journal of Universal Computer Science
}
\author{Andrius Merkys}
