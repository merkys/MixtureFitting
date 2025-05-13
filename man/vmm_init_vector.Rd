\name{vmm_init_vector}
\alias{vmm_init_vector}
\docType{methods}
\title{Estimate von Mises Mixture parameters using Expectation Maximisation.}
\description{
    Estimate an initialization vector for von Mises mixture fitting via Expectation Maximisation.
    Proportions are set to equal, centers are equispaced through the whole domain of input sample, and concentrations are set to (m/(12*180))^2.
}
\usage{
    vmm_init_vector( m, implementation = "C" )
}
\arguments{
    \item{m}{number of mixture components}
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{
    Parameter vector of 3*\emph{n} parameters, where \emph{n} is
    number of mixture components. Structure of p vector is
    p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n} ),
    where A\emph{i} is the proportion of \emph{i}-th component,
    mu\emph{i} is the center of \emph{i}-th component and
    k\emph{i} is the concentration of \emph{i}-th component.
}
\author{Andrius Merkys}
