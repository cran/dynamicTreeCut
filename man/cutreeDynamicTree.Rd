\name{cutreeDynamicTree}
\alias{cutreeDynamicTree}
\title{Dynamic Dendrogram Pruning Based on Dendrogram Only}
\description{
  Detect clusters in a hierarchical dendrogram using a variable cut height approach. Uses only the
information in the dendrogram itself is used (which may give incorrect assignment for outlying objects).
}
\usage{
cutreeDynamicTree(dendro, maxTreeHeight = 1, deepSplit = TRUE, minModuleSize = 50)
}
\arguments{
  \item{dendro}{Hierarchical clustering dendrogram such produced by \code{\link[stats]{hclust}}.}
  \item{maxTreeHeight}{Maximum joining height of objects to be considered part of clusters.}
  \item{deepSplit}{If \code{TRUE}, method will favor sensitivity and produce more smaller clusters. When
\code{FALSE}, there will be fewer bigger clusters.}
  \item{minModuleSize}{Minimum module size. Branches containing fewer than \code{minModuleSize} objects
will be left unlabeled.}
}
\details{
  A variable height branch pruning technique for dendrograms produced by hierarchical clustering.
Initially, branches are cut off at the height \code{maxTreeHeight}; the resulting clusters are then
examined for substructure and if subclusters are detected, they are assigned separate labels. Subclusters
are detected by structure and are required to have a minimum of \code{minModuleSize} objects on them to
be assigned a separate label. A rough degree of control over what it means to be a subcluster is
implemented by the parameter \code{deepSplit}.
}
\value{
  A vector of numerical labels giving assignment of objects to modules. Unassigned objects are labeled 0,
the largest module has label 1, next largest 2 etc.
}
\references{ 
  \url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting}
}
\author{ Bin Zhang, \email{binzhang.ucla@gmail.com}, with contributions by Peter Langfelder,
\email{Peter.Langfelder@gmail.com}.}
\seealso{\code{\link[stats]{hclust}}, \code{\link{cutreeHybrid}}}
\keyword{misc}
