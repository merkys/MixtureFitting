\name{smm_fit_em_CWL04}
\alias{smm_fit_em_CWL04}
\docType{methods}
\title{Greedily estimate Student's t Mixture parameters using Expectation Maximization.}
\description{
    Estimates (greedily) parameters for univariate Student's t mixture
    using Expectation Maximization algorithm, implemented according to
    Chen et al. (2004). The algorithm relies upon smm_fit_em_GNL08()
    to estimate mixture parameters iteratively.
}
\usage{
    smm_fit_em_CWL04( x, p, collect.history = FALSE, debug = FALSE,
                      ... )
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
    \item{collect.history}{
        logical. If set to TRUE, a list of parameter values of all
        iterations is returned.
    }
    \item{debug}{
        flag to turn the debug prints on/off.
    }
    \item{...}{
        parameters passed to smm_fit_em_GNL08().
    }
}
\value{
    A list.
}
\references{
    Chen, S.; Wang, H. & Luo, B.
    Greedy EM Algorithm for Robust T-Mixture Modeling
    Third International Conference on Image and Graphics (ICIG'04),
    Institute of Electrical & Electronics Engineers (IEEE), 2004, 548--551
}
\author{Andrius Merkys}
