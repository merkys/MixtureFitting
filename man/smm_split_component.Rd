\name{smm_split_component}
\alias{smm_split_component}
\docType{methods}
\title{Split a component of Student's t-distribution in two.}
\description{
    Splits a component of Student's t-distribution mixture. Implemented
    according to Eqns. 30--36 of Chen et al. (2004).
}
\usage{
    smm_split_component( p, alpha = 0.5, beta = 0.5, u = 0.5 )
}
\arguments{
    \item{p}{
        vector of Student's t mixture parameters. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n}, ni1, ni2, ..., ni\emph{n} ),
        where \emph{n} is number of mixture components,
        A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        k\emph{i} is the concentration of \emph{i}-th component and
        ni\emph{i} is the degrees of freedom of \emph{i}-th component.
    }
    \item{alpha}{split proportion for component proportions}
    \item{beta}{split proportion for component concentrations}
    \item{u}{split proportion for component centers}
}
\value{
    Vector of parameters for resulting two-component mixture, whose
    structure is the same as of input parameter's p.
}
\references{
    Chen, S.-B. & Luo, B.
    Robust t-mixture modelling with SMEM algorithm
    Proceedings of 2004 International Conference on Machine Learning and
    Cybernetics (IEEE Cat. No.04EX826),
    Institute of Electrical & Electronics Engineers (IEEE), 2004, 6, 3689--3694
}
\author{Andrius Merkys}
