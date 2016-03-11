\name{merge2Clusters}
\Rdversion{1.1}
\alias{merge2Clusters}
\title{
  Merge Two Clusters
}
\description{
 Merge 2 clusters into 1.
}
\usage{
merge2Clusters(labels, mainClusterLabel, minorClusterLabel)
}
\arguments{

  \item{labels}{
  a vector or factor giving the cluster labels
}
  \item{mainClusterLabel}{
  label of the first merged cluster. The merged cluster will have this label.
}
  \item{minorClusterLabel}{
  label of the second merged cluster.
}
}
\value{
  A vector or factor of the merged labels.
}
\author{
 Bin Zhang and Peter Langfelder
}
\examples{

options(stringsAsFactors = FALSE);

# Works with character labels:
labels = c(rep("grey", 5), rep("blue", 2), rep("red", 3))
merge2Clusters(labels, "blue", "red")

# Works with factor labels:
labelsF = factor(labels)
merge2Clusters(labelsF, "blue", "red")

# Works also with numeric labels:

labelsN = as.numeric(factor(labels))
labelsN
merge2Clusters(labelsF, 1, 3)


}
\keyword{ misc }
