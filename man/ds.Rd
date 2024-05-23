\name{ds}
\alias{ds}
\docType{methods}
\title{Density of The Student's t Model}
\description{
    Density function for the Student's t Model.
    Wrapper around R's dt(), supporting center and concentration parameters.
}
\usage{
    ds( x, c, s, ni )
}
\arguments{
    \item{x}{data vector}
    \item{c}{center}
    \item{s}{concentration}
    \item{ni}{degrees of freedom}
}
\value{A vector.}
\author{Andrius Merkys}
