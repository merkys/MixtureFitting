\name{gmm_merge_components}
\alias{gmm_merge_components}
\docType{methods}
\title{Merge two Gaussian components into one.}
\description{
    Merges \emph{i}th and \emph{j}th components of Gaussian mixture
    model. Implemented in the same venue as in
    \code{\link[fpc]{mergeparameters}} of \code{fpc}.
}
\usage{
    gmm_merge_components( x, p, i, j )
}
\arguments{
    \item{x}{data vector}
    \item{p}{
        vector of Gaussian mixture parameters. Structure of p vector is
        p = c( A1, A2, ..., A\emph{n}, mu1, mu2, ..., mu\emph{n}, k1, k2, ..., k\emph{n} ),
        where \emph{n} is number of mixture components,
        A\emph{i} is the proportion of \emph{i}-th component,
        mu\emph{i} is the center of \emph{i}-th component,
        k\emph{i} is the concentration of \emph{i}-th component.
    }
    \item{i}{
        index of the first component to be merged. Component with this
        index will be replaced by a merged one in the output.
    }
    \item{j}{
        index of the second component to be merged. Component with this
        index will be removed in the output.
    }
}
\value{
    Vector of mixture parameters, whose structure is the same as of input
    parameter's p.
}
\references{
    Hennig, C.
    Methods for merging Gaussian mixture components
    Advances in Data Analysis and Classification, Springer Nature, 2010, 4, 3-34
}
\author{Andrius Merkys}
