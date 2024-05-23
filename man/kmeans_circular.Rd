\name{kmeans_circular}
\alias{kmeans_circular}
\docType{methods}
\title{K-Means Clustering for Points on Circle}
\description{
    Perform k-means clustering on angular data (in degrees).
}
\usage{
    kmeans_circular( x, centers, iter.max = 10 )
}
\arguments{
    \item{x}{data vector}
    \item{centers}{vector of initial centers (in degrees)}
    \item{iter.max}{maximum number of iterations}
}
\value{
    Vector of the same length as \emph{centers} defining cluster centers (in degrees).
}
\author{Andrius Merkys}
