\name{smm_init_vector}
\alias{smm_init_vector}
\docType{methods}
\title{Estimate Student's t Mixture parameters using Expectation Maximization.}
\description{
    Estimate an initialization vector for Student's t mixture fitting via Expectation Maximization.
    Proportions are set to be equal, centers are equispaced through the whole domain of input sample, concentrations and degrees of freedom are set to 1.
}
\usage{
    smm_init_vector( x, n )
}
\arguments{
    \item{x}{data vector}
    \item{n}{number of mixture components}
}
\value{
    Parameter vector of 4*\emph{n} parameters, where \emph{n} is number of mixture components.
    Structure of p vector is
    p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n}, ni1, ni2, ..., ni\emph{n} ),
    where A\emph{i} is the proportion of \emph{i}-th component,
    mu\emph{i} is the center of \emph{i}-th component,
    k\emph{i} is the concentration of \emph{i}-th component and
    ni\emph{i} is the degrees of freedom of \emph{i}-th component.
}
\author{Andrius Merkys}
