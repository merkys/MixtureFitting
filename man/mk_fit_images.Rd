\name{mk_fit_images}
\alias{mk_fit_images}
\docType{methods}
\title{Mixture Distribution Modelling}
\description{
    Draw a PNG histogram with a mixture density on top of it for each iteration of mixture optimization process.
}
\usage{
    mk_fit_images( h, l, prefix = "img_" )
}
\arguments{
    \item{h}{histogram object, as returned from hist()}
    \item{l}{list containing model vectors}
    \item{prefix}{prefix of file name to write}
}
\author{Andrius Merkys}
