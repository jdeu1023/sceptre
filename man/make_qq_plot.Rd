% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/high_moi.R
\name{make_qq_plot}
\alias{make_qq_plot}
\title{Make QQ-plot}
\usage{
make_qq_plot(p_values, ci_level = 0.95, point_col = "royalblue4", alpha = 0.9)
}
\arguments{
\item{p_values}{a vector values -- hypothesized to follow a uniform distribution -- from which to construct the QQ-plot. This vector typically will be a set of \emph{p}-values.}

\item{ci_level}{level of the pointwise confidence band (default 0.95)}

\item{point_col}{color of the plotted points}

\item{alpha}{transparency of the plotted points}
}
\value{
a ggplot object containing the QQ-plot
}
\description{
Makes a QQ-plot for a set of values hypothesized to follow a uniform distribution (e.g., \emph{p}-values).
}
\examples{
set.seed(4)
p_vals <- runif(5000)
make_qq_plot(p_vals)
}
