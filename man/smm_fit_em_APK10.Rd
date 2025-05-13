\name{smm_fit_em_APK10}
\alias{smm_fit_em_APK10}
\docType{methods}
\title{Estimate Student's t Mixture parameters using Expectation Maximisation.}
\description{
    Estimates parameters for univariate Student's t mixture using
    Expectation Maximisation algorithm, according to Fig. 2 of
    Aeschliman et al. (2010).
}
\usage{
    smm_fit_em_APK10( x, p, epsilon = c( 1e-6, 1e-6, 1e-6, 1e-6 ),
                      collect.history = FALSE, debug = FALSE )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        initialization vector of 4*\emph{n} parameters, where \emph{n} is
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
    \item{collect.history}{
        flag to turn accumulation of estimation history on/off.
    }
    \item{debug}{
        flag to turn the debug prints on/off.
    }
}
\value{
    A list.
}
\references{
    Aeschliman, C.; Park, J. & Kak, A. C. A
    Novel Parameter Estimation Algorithm for the Multivariate t-Distribution
    and Its Application to Computer Vision
    European Conference on Computer Vision 2010, 2010
    \url{https://engineering.purdue.edu/RVL/Publications/Aeschliman2010ANovel.pdf}
}
\author{Andrius Merkys}
