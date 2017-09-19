\name{kldiv}
\alias{kldiv}
\docType{methods}
\title{Kullback--Leibler Divergence of \emph{i}th Student's t Mixture component.}
\description{
    Measures Kullback--Leibler divergence of \emph{i}th Student's t
    Mixture component using Dirac's delta function. Implemented
    according to Chen et al. (2004).
}
\usage{
    kldiv( x, p, k )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        vector of Student's t mixture parameters. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n}, ni1, ni2, ..., ni\emph{n} ),
        where \emph{n} is number of mixture components,
        A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        k\emph{i} is the concentration of \emph{i}-th component and
        ni\emph{i} is the degrees of freedom of \emph{i}-th component.
    }
    \item{k}{
        number of the component.
    }
}
\value{
    Kullback--Leibler divergence as double.
}
\references{
    Chen, S.; Wang, H. & Luo, B.
    Greedy EM Algorithm for Robust T-Mixture Modeling
    Third International Conference on Image and Graphics (ICIG'04),
    Institute of Electrical & Electronics Engineers (IEEE), 2004, 548--551
}
\author{Andrius Merkys}
