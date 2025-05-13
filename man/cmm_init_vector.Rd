\name{cmm_init_vector}
\alias{cmm_init_vector}
\docType{methods}
\title{Estimate Cauchy Mixture parameters using Expectation Maximisation.}
\description{
    Estimate an initialization vector for Cauchy mixture fitting via Expectation Maximisation.
    Proportions are set to equal, centers are equispaced through the whole domain of input sample, and scales are set to 1.
}
\usage{
    cmm_init_vector( x, m, implementation = "C" )
}
\arguments{
    \item{x}{data vector}
    \item{m}{number of mixture components}
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
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
