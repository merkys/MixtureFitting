\name{cmm_fit_em}
\alias{cmm_fit_em}
\docType{methods}
\title{Estimate Cauchy Mixture parameters using Expectation Maximisation.}
\description{
    Estimates parameters for Caucy mixture using Expectation Maximisation algorithm.
}
\usage{
    cmm_fit_em( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ),
                iter.cauchy = 20, debug = FALSE, implementation = "C" )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialization vector of 3*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, gamma1, gamma2, ..., gamma\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component and
        gamma\emph{i} is the Cauchy scale of \emph{i}-th component.
    }
    \item{epsilon}{
        tolerance threshold for convergence. Structure of epsilon is
        epsilon = c( epsilon_A, epsilon_mu, epsilon_gamma ), where
        epsilon_A is threshold for component proportions,
        epsilon_mu is threshold for component centers and
        epsilon_gamma is threshold for component Cauchy scales.
    }
    \item{iter.cauchy}{
        number of iterations to fit a single Cauchy component.
    }
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
    Ferenc Nahy. Parameter Estimation of the Cauchy Distribution in
    Information Theory Approach (2006).
    Journal of Universal Computer Science
}
\author{Andrius Merkys}
