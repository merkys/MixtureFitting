\name{smm_fit_em_GNL08}
\alias{smm_fit_em_GNL08}
\docType{methods}
\title{Estimate Student's t Mixture parameters using Expectation Maximisation.}
\description{
    Estimates parameters for univariate Student's t mixture using
    Expectation Maximisation algorithm, according to Eqns. 12--17 of
    Gerogiannis et al. (2009).
}
\usage{
    smm_fit_em_GNL08( x, p, epsilon = c( 1e-6, 1e-6, 1e-6, 1e-6 ),
                      collect.history = FALSE, debug = FALSE,
                      min.sigma = 1e-256, min.ni = 1e-256,
                      max.df = 1000, max.steps = Inf,
                      polyroot.solution = 'jenkins_taub',
                      convergence = abs_convergence,
                      unif.component = FALSE )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialisation vector of 4*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n}, ni1, ni2, ..., ni\emph{n} ),
        where A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        k\emph{i} is the concentration of \emph{i}-th component and
        ni\emph{i} is the degrees of freedom of \emph{i}-th component.
    }
    \item{epsilon}{
        tolerance threshold for convergence. Structure of epsilon is
        epsilon = c( epsilon_A, epsilon_mu, epsilon_k, epsilon_ni ), where
        epsilon_A is threshold for component proportions,
        epsilon_mu is threshold for component centers,
        epsilon_k is threshold for component concentrations and
        epsilon_ni is threshold for component degrees of freedom.
    }
    \item{debug}{
        flag to turn the debug prints on/off.
    }
    \item{max.steps}{maximum number of steps, may be infinity}
    \item{convergence}{
        function to use for convergence checking.
        Must accept function values of the last two iterations and return TRUE or FALSE.
    }
    \item{unif.component}{should a uniform component for outliers be added?}
}
\value{
    A list.
}
\references{
    Gerogiannis, D.; Nikou, C. & Likas, A.
    The mixtures of Student's t-distributions as a robust framework
    for rigid registration.
    Image and Vision Computing, Elsevier BV, 2009, 27, 1285--1294
    \url{http://www.cs.uoi.gr/~arly/papers/imavis09.pdf}
}
\author{Andrius Merkys}
