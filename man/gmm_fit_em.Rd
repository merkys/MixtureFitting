\name{gmm_fit_em}
\alias{gmm_fit_em}
\docType{methods}
\title{Estimate Gaussian Mixture parameters using Expectation Maximization.}
\description{
    Estimates parameters for Gaussian mixture using Expectation Maximization algorithm.
}
\usage{
    gmm_fit_em( x, p, w = numeric(), epsilon = c( 0.000001, 0.000001, 0.000001 ),
                debug = FALSE, implementation = "C", ... )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialization vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, sigma1, sigma2, ..., sigma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        sigma\emph{i} is the scale of \emph{i}-th component.
    }
    \item{w}{
        weights of data points, must have the same length as the data
        vector; if not given or has different length, equal weights are
        assumed.
    }
    \item{epsilon}{
        tolerance threshold for convergence. Structure of epsilon is
        epsilon = c( epsilon_A, epsilon_mu, epsilon_sigma ), where
        epsilon_A is threshold for component proportions,
        epsilon_mu is threshold for component centers and
        epsilon_sigma is threshold for component scales.
    }
    \item{debug}{
        flag to turn the debug prints on/off.
    }
    \item{implementation}{
        flag to switch between C (default) and R implementations.
    }
    \item{...}{
        additional arguments passed to gmm_fit_em_R() when R implementation
        is used.
    }
}
\value{
    Vector of mixture parameters, whose structure is the same as of input
    parameter's p.
}
\author{Andrius Merkys}
