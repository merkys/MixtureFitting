\name{rcmm}
\alias{rcmm}
\docType{methods}
\title{Random Sample of The Cauchy Mixture Distribution}
\description{
    Generates a random sample of the Cauchy mixture distribution.
}
\usage{
    rcmm( n, p )
}
\arguments{
    \item{n}{sample size}
    \item{p}{
        parameter vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, gamma1, gamma2, ..., gamma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the location of \emph{i}-th component,
        gamma\emph{i} is the Cauchy scale of \emph{i}-th component.
    }
}
\value{A vector.}
\author{Andrius Merkys}
