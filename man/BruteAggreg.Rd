\name{BruteAggreg}
\alias{BruteAggreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Rank Aggregation via brute force algorithm}
\description{
 Weighted rank aggregation of ordered lists is performed using the brute force approach, i.e. 
 generating all possible ordered lists and finding the list with the minimum value of the objective
 function
}
\usage{
BruteAggreg(x, k, weights = NULL, distance = c("Spearman", "Kendall"),
            importance=rep(1,nrow(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix of ordered lists to be combined (lists must be in rows)}
  \item{k}{size of the top-k list}
  \item{weights}{scores (weights) to be used in the aggregation process}
  \item{distance}{distance which "measures" the similarity between the ordered lists}
  \item{importance}{a vector of weights indicating the importance of each ordered list in x}
}
\details{
  The function performs rank aggregation using the old-fashion brute force approach. This approach
  works for small problems only and should not be attempted if k is relatively large (k > 10). To 
  generate all possible ordered lists, the permutation function from the gtools package is used. Both
  weighted and unweighted rank aggregation can be performed. Please refer to the documentation for 
  RankAggreg function as the same constraints on x and index.weights apply to both functions.
}
\value{
  \item{top.list}{Top-k aggregated list}
  \item{optimal.value}{the minimum value of the objective function corresponding to the top-k list}
  \item{distance}{distance used by the algorithm}
  \item{method}{method used: BruteForce}
  \item{importance}{importance vector used}
  \item{lists}{original lists to be combined}
  \item{weights}{scaled weights used in aggregation}
  \item{sample}{objective function values}
  \item{sample.size}{number of all possible solutions}
  \item{summary}{contains minimum and median values of sample}
}

\references{Pihur, V., Datta, S., and Datta, S. (2007) "Weighted rank aggregation of cluster validation 
measures: a Monte Carlo cross-entropy approach" Bioinformatics, 23(13):1607-1615 }

\author{Vasyl Pihur, Somnath Datta, Susmita Datta}

\seealso{\code{\link{RankAggreg}}}

\examples{
require(gtools)

# rank aggregation without weights
x <- matrix(c("A", "B", "C", "D", "E",
        "B", "D", "A", "E", "C",
        "B", "A", "E", "C", "D",
        "A", "D", "B", "C", "E"), byrow=TRUE, ncol=5)

(toplist <- BruteAggreg(x, 5))

# weighted rank aggregation
set.seed(100)
w <- matrix(rnorm(20), ncol=5)
w <- t(apply(w, 1, sort))

(toplist <- BruteAggreg(x,5,w,"Spearman")) # using the Spearman distance
(toplist <- BruteAggreg(x,5,w,"Kendall")) #using the Kendall distance
plot(toplist)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}
\keyword{robust}% __ONLY ONE__ keyword per line
