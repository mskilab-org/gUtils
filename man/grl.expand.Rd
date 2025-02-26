% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gUtils.R
\name{grl.expand}
\alias{grl.expand}
\title{grl.expand}
\usage{
grl.expand(grl, expand_win)
}
\arguments{
\item{grl}{\code{GRangesList}}

\item{expand_win}{\code{integer}}
}
\value{
GRangesList with added window
}
\description{
Function wrapping around the `+` operator
for GRanges objects to work on GRangesLists.
Expands window of element GRanges within GRangesList
}
\author{
Kevin Hadi
}
