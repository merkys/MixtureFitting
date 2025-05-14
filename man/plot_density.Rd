\name{plot_density}
\alias{plot_density}
\docType{methods}
\title{Mixture Distribution Modeling}
\description{
    Draw a PNG histogram with a mixture density on top of it.
}
\usage{
    plot_density( x, model, density_f, width, height,
                  cuts = 400, main = '',
                  filename = '/dev/stdout',
                  obs_good = c(), obs_bad = c(),
                  scale_density = FALSE )

}
\arguments{
    \item{x}{data vector}
    \item{cuts}{number of breaks in histogram}
    \item{main}{main title of the plot}
    \item{model}{model passed to density_f()}
    \item{density_f}{probability density function}
    \item{filename}{name of the file to write}
    \item{width}{image width, passed to png()}
    \item{height}{image height, passed to png()}
    \item{obs_good}{vector of values to mark with rug() in green color}
    \item{obs_bad}{vector of values to mark with rug() in red color}
    \item{scale_density}{should probability density be scaled?}
}
\author{Andrius Merkys}
