\name{vmm_fit_em_by_ll}
\alias{vmm_fit_em_by_ll}
\docType{methods}
\title{Estimate von Mises Mixture parameters using Expectation Maximisation.}
\description{
    Estimates parameters for univariate von Mises mixture using Expectation Maximisation algorithm.
    In this version stopping criterion is the difference between log-likelihood estimates of subsequent iterations.
}
\usage{
    vmm_fit_em_by_ll( x, p, epsilon = .Machine$double.eps,
                      debug = FALSE, implementation = "C" )
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
    \item{epsilon}{tolerance threshold for convergence}
    \item{debug}{
        flag to turn the debug prints on/off.
    }
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
}
\value{
    Vector of mixture parameters, whose structure is the same as of input
    parameter's p.
}
\references{
    Banerjee et al.
    Expectation Maximization for Clustering on Hyperspheres (2003),
    manuscript, accessible on:
    \url{https://web.archive.org/web/20130120061240/http://www.lans.ece.utexas.edu/~abanerjee/papers/05/banerjee05a.pdf}
}
\author{Andrius Merkys}
