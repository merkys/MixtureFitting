\name{snmm_fit_em}
\alias{snmm_fit_em}
\docType{methods}
\title{Estimate Skew-Normal Mixture parameters using Expectation Maximization.}
\description{
    Estimates parameters for skew-normal mixture using Expectation Maximization algorithm.
}
\usage{
    snmm_fit_em( x, p, epsilon = c( 0.000001, 0.000001, 0.000001 ) )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialization vector of 4*\emph{n} parameters, where \emph{n} is
        number of mixture components. Structure of p vector is
        p = c( omega1, omega2, ..., omega\emph{n}, dzeta1, dzeta2, ..., dzeta\emph{n}, sigma1, sigma2, ..., sigma\emph{n}, lambda1, lambda2, ..., lambda\emph{n} ),
        where omega\emph{i} is the proportion of \emph{i}-th component,
        dzeta\emph{i} is the location of \emph{i}-th component,
        sigma\emph{i} is the scale of \emph{i}-th component and
        lambda\emph{i} is the shape of \emph{i}-th component.
    }
    \item{epsilon}{
        tolerance threshold for convergence. Structure of epsilon is
        epsilon = c( epsilon_omega, epsilon_dzeta, epsilon_sigma, epsilon_lambda ), where
        epsilon_omega is threshold for component proportions,
        epsilon_dzeta is threshold for component locations,
        epsilon_sigma is threshold for component scales and
        epsilon_lambda is threshold for component shapes.
    }

}
\value{
    Vector of mixture parameters, whose structure is the same as of input
    parameter's p.
}
\references{
    Arnold et al.
    The nontruncated marginal of a truncated bivariate normal distribution,
    Psychometrika 58, pages 471--488 (1993)

    Lin et al.
    Finite mixture modelling using the skew normal distribution,
    Statistica Sinica 17 (2007), 909--927
    \url{https://www3.stat.sinica.edu.tw/statistica/oldpdf/A17n35.pdf}
}
\author{Andrius Merkys}
