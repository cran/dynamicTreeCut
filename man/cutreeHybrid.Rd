\name{cutreeHybrid}
\alias{cutreeHybrid}
\title{Hybrid adaptive tree cut for hierarchical clustering dendrograms.  }
\description{
  Detect clusters in a dendorgram produced by the function \code{hclust}.  
}
\usage{
cutreeHybrid(dendro, distM, cutHeight = NULL, minClusterSize = 20, 
             deepSplit = 1,
             maxCoreScatter = NULL, minGap = NULL,
             maxAbsCoreScatter = NULL, minAbsGap = NULL, 
             pamStage = TRUE, pamRespectsDendro = TRUE,
             useMedoids = FALSE, maxDistToLabel = cutHeight,
             respectSmallClusters = TRUE, verbose = 2, indent = 0)
}
\arguments{
  \item{dendro}{a hierarchical clustering dendorgram such as one returned by \code{hclust}.  }
  \item{distM}{Distance matrix that was used as input to \code{hclust}.  }
  \item{cutHeight}{Maximum joining heights that will be considered. It defaults to 99% 
of the range between the 5th percentile and the maximum of the joining heights on the dendrogram.  }
  \item{minClusterSize}{Minimum cluster size.  }
  \item{deepSplit}{Either logical or integer in the range 0 to 4. Provides a rough control over 
sensitivity to cluster splitting. The higher the value, the more and smaller clusters will be produced. 
A finer control can be achieved via \code{maxCoreScatter} and \code{minGap} below.}
  \item{maxCoreScatter}{ Maximum scatter of the core for a branch to be a cluster, given as the fraction
of \code{cutHeight} relative to the 5th percentile of joining heights. See Details.  } 
  \item{minGap}{Minimum cluster gap given as the fraction of the difference between \code{cutHeight} and
the 5th percentile of joining heights.  }
  \item{maxAbsCoreScatter}{Maximum scatter of the core for a branch to be a cluster given as absolute
heights. If given, overrides \code{maxCoreScatter}.  }
  \item{minAbsGap}{ Minimum cluster gap given as absolute height difference. If given, overrides
\code{minGap}.  } 
  \item{pamStage}{Logical, only used for method "hybrid".  If \code{TRUE}, the second (PAM-like) stage will be performed.
}    
  \item{pamRespectsDendro}{Logical, only used for method "hybrid". 
  If \code{TRUE}, the PAM stage will
respect the dendrogram in the sense an object can be PAM-assigned only to clusters that lie below it on 
the branch that the object is merged into.
     See \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details. }

\item{useMedoids}{if TRUE, the second stage will be use object to medoid distance; if FALSE, it
  will use average object to cluster distance. The default (FALSE) is recommended. }
  \item{maxDistToLabel}{Maximum object distance to closest cluster that will result in the object
assigned to that cluster.  }
  \item{respectSmallClusters}{If TRUE, branches that failed to be clusters in stage 1 only because of
insufficient size will be assigned together in stage 2. If FALSE, all objects will be assigned
individually.}
  \item{verbose}{Controls the verbosity of the output. 0 will make the function completely quiet,
values up to 4 gradually increase verbosity.}
  \item{indent}{Controls indentation of printed messages (see \code{verbose} above). Each unit
adds two spaces before printed messages; useful when several functions' output is to be nested. } 
}
\details{
  The function detects clusters in a hierarchical dendrogram based on the shape of branches on the
dendrogram. For details on the method, see
\url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting}.

  In order to make the shape parameters \code{maxCoreScatter} and \code{minGap} more universal, their
values are interpreted relative to \code{cutHeight} and the 5th percetile of the merging heights (we
arbitrarily chose the 5th percetile rather than the minimum for reasons of stability). Thus, the absolute
maximum allowable core scatter is calculated as \code{maxCoreScatter * (cutHeight - refHeight) +
refHeight} and the absolute minimum allowable gap as \code{minGap * (cutHeight - refHeight)}, where
\code{refHeight} is the 5th percentile of the merging heights.

}
\value{
  A list containg the following elements:
  \item{labels}{Numerical labels of clusters, with 0 meaning unassigned, label 1 labeling the largest
                cluster etc.  }
  \item{cores}{Numerical labels indicating cores of found clusters.  }
  \item{smallLabels}{Numerical labels for branches that failed to be recognized clusters only
because of insufficient number of objects.}
  \item{trimmed}{Numerical labels indicating objects that have been trimmed from clusters.  }
  \item{branches}{A list detailing the deteced branch structure.  }
}
\references{
  Langfelder P, Zhang B, Horvath S, 2007.
  \url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting}
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{
  \code{\link[stats]{hclust}}, \code{\link[stats]{as.dist}} 
}
\keyword{cluster}
