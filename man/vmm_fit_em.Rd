\name{vmm_fit_em}
\alias{vmm_fit_em}
\docType{methods}
\title{Estimate von Mises Mixture parameters using Expectation Maximisation.}
\description{
    Estimates parameters for univariate von Mises mixture using Expectation
    Maximisation algorithm.
}
\usage{
    vmm_fit_em( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ) )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialisation vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        k\emph{i} is the concentration of \emph{i}-th component.
    }
    \item{epsilon}{
        tolerance threshold for convergence. Structure of epsilon is
        epsilon = c( epsilon_A, epsilon_mu, epsilon_k ), where
        epsilon_A is threshold for component proportions,
        epsilon_mu is threshold for component centers and
        epsilon_k is threshold for component concentrations.
    }
}
\value{
    Vector of mixture parameters, whose structure is the same as of input
    parameter's p.
}
\references{
    Banerjee et al.
    Expectation Maximization for Clustering on Hyperspheres (2003),
    manuscript, accesible on:
    \url{http://www.lans.ece.utexas.edu/~abanerjee/papers/05/banerjee05a.pdf}
}
\author{Andrius Merkys}
