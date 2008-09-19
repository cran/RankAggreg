\name{geneLists}
\alias{geneLists}
\docType{data}
\title{Ordered Gene Lists from 5 microarray studies}
\description{
  This dataset contains five lists of genes, each of size 25, from five independent
  microarray studies on prostate cancer. The lists are given in Table 4 in the 
  manuscript by DeConde et al. Lists form the rows of the dataset with columns corresponding to the ranks of 
   genes in each individual study.
}
\usage{data(geneLists)}
\format{
  A matrix of size 5 by 25 containing 5 lists of genes.
}
\source{
  R. DeConde, S. Hawley, S. Falcon, N. Clegg, B. Knudsen, and R. Etzioni. 
Combining results of microarray experiments: a rank aggregation approach. Stat Appl
Genet Mol Biol, 5(1):Article15, 2006.
}

\examples{
data(geneLists)
topList <- RankAggreg(geneLists, 5, N=700, seed=100, convIn=3)
plot(topList)
}
\keyword{datasets}
