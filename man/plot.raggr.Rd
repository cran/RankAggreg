\name{plot.raggr}
\alias{plot.raggr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot function for raggr object returned by RankAggreg or BruteAggreg}
\description{
  Plots individual ordered lists with the corresponding solution. Optionally,
  naive average rank aggregation can be added.
}
\usage{
plot.raggr(x, show.average = TRUE, show.legend = TRUE, colR="red", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{raggr object returned by RankAggreg}
  \item{show.average}{boolean if average aggregation to be plotted}
  \item{show.legend}{boolean if the legend is to be displayed}
  \item{colR}{specifies the color for the resulting list}
  \item{\dots}{additional plotting parameters}
}
\details{
  The function plots individual lists and the solution using ranks only (weights are
  not used at any time). Optional average rank aggregation can be performed and visualized.
  Average rank aggregation is a simple aggregation procedure which computes the average ranks
  for each unique element accross and orders them from the smallest to the largest value.
}
\value{
   Nothing is returned
}
\references{Pihur, V., Datta, S., and Datta, S. (2007) "Weighted rank aggregation of cluster validation 
measures: a Monte Carlo cross-entropy approach" Bioinformatics, 23(13):1607-1615 }

\author{Vasyl Pihur, Somnath Datta, Susmita Datta}

\seealso{\code{\link{RankAggreg}}, \code{\link{BruteAggreg}}}
\examples{
# rank aggregation without weights
x <- matrix(c("A", "B", "C", "D", "E",
        "B", "D", "A", "E", "C",
        "B", "A", "E", "C", "D",
        "A", "D", "B", "C", "E"), byrow=TRUE, ncol=5)

(CES <- RankAggreg(x, 5, method="CE", distance="Spearman", rho=.1, verbose=FALSE))
plot(CES)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
\keyword{robust}% __ONLY ONE__ keyword per line
