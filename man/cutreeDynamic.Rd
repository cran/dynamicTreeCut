\name{cutreeDynamic}
\alias{cutreeDynamic}
\title{Adaptive branch pruning of hierarchical clustering dendrograms.}
\description{
  This wrapper provides a common access point for two methods of adaptive branch pruning of hierarchical
clustering dendrograms. 
}
\usage{
cutreeDynamic(
      dendro, cutHeight = NULL, minClusterSize = 20,

      # Basic tree cut options
      method = "hybrid", distM = NULL,
      deepSplit = (ifelse(method=="hybrid", 1, FALSE)),

      # Advanced options
      maxCoreScatter = NULL, minGap = NULL,
      maxAbsCoreScatter = NULL, minAbsGap = NULL,

      # External (user-supplied) measure of branch split
      externalBranchSplitFnc = NULL, minExternalSplit = NULL,
      externalSplitOptions = list(),
      assumeSimpleExternalSpecification = TRUE,

      # PAM stage options
      pamStage = TRUE, pamRespectsDendro = TRUE,
      useMedoids = FALSE, maxDistToLabel = NULL,
      maxPamDist = cutHeight,
      respectSmallClusters = TRUE,

      # Various options
      verbose = 2, indent = 0)

}
\arguments{
  \item{dendro}{A hierarchical clustering dendorgram such as one returned by \code{hclust}.  }

  \item{cutHeight}{Maximum joining heights that will be considered. For \code{method=="tree"} it
defaults to 0.99. For \code{method=="hybrid"} it defaults to 99\% of the range between the 5th
percentile and the maximum of the joining heights on the dendrogram.}

  \item{minClusterSize}{Minimum cluster size.  }

  \item{method}{Chooses the method to use. Recognized values are "hybrid" and "tree".  }

  \item{distM}{Only used for method "hybrid". The distance matrix used as input to \code{hclust}. If not
given and \code{method == "hybrid"}, the function will issue a warning and default to \code{method =
"tree"}.}

  \item{deepSplit}{For method "hybrid", can be either logical or integer in the range 0 to 4. For method
"tree", must be logical. In both cases, provides a rough control over sensitivity to cluster splitting.
The higher the value (or if \code{TRUE}), the more and smaller clusters will be produced. For the
"hybrid" method, a finer control can be achieved via \code{maxCoreScatter} and \code{minGap} below.}

  \item{maxCoreScatter}{Only used for method "hybrid". 
Maximum scatter of the core for a branch to be a cluster, given as the fraction of \code{cutHeight}
relative to the 5th percentile of joining heights. See Details.  }

  \item{minGap}{Only used for method "hybrid". 
Minimum cluster gap given as the fraction of the difference between \code{cutHeight} and the 5th
percentile of joining heights.  }

  \item{maxAbsCoreScatter}{Only used for method "hybrid". 
Maximum scatter of the core for a branch to be a cluster given as absolute heights. If given, overrides
\code{maxCoreScatter}.  }

  \item{minAbsGap}{Only used for method "hybrid". 
Minimum cluster gap given as absolute height difference. If given, overrides \code{minGap}.  }

  \item{externalBranchSplitFnc}{Optional function to evaluate split (dissimilarity) between two branches.
Either a single function or a list in which each component is a function (see
\code{assumeSimpleExternalSpecification} below for how to specify a single function).
Each function can be specified by name (a character string) or the actual function object.
Each given function must take arguments \code{branch1} and \code{branch2} that specify the 
indices of objects in
the two branches whose dissimilarity is to be evaluated, and possibly other arguments. It must return a
number that quantifies the dissimilarity of the two branches. The returned value will be compared to
\code{minExternalSplit} (see below). This argument is only used for  method "hybrid".}

  \item{minExternalSplit}{Thresholds to decide whether two branches should be merged. 
It should be a numeric vector of the same length as the number of functions in
\code{externalBranchSplitFnc} above. 
Only used for method "hybrid". }

  \item{externalSplitOptions}{Further arguments to function \code{externalBranchSplitFnc}. If only one
external function is specified in \code{externalBranchSplitFnc} above, \code{externalSplitOptions} 
can be a named
list of arguments or a list with one component that is in turn  the named list of further arguments for
\code{externalBranchSplitFnc[[1]]}. The argument \code{assumeSimpleExternalSpecification} controls which of
the two possibilities should be assumed. If multiple functions are specified in \code{externalBranchSplitFnc},
\code{externalSplitOptions}
must be a list in which each component is a named list giving the further arguments for the corresponding 
function in  \code{externalBranchSplitFnc}. Only used for method "hybrid".}

  \item{assumeSimpleExternalSpecification}{Logical: when \code{minExternalSplit} above is a scalar (has
length 1), should the function assume a simple specification of \code{externalBranchSplitFnc} and
\code{externalSplitOptions}? If \code{TRUE}, \code{externalBranchSplitFnc} is taken as the function
specification and \code{externalSplitOptions} the named list of options. This is suitable for simple direct
calls of this function. If \code{FALSE}, \code{externalBranchSplitFnc} is assumed to be a list with a
single component which specifies the function, and \code{externalSplitOptions} is a list with one
component that is in turn the named list of further arguments for \code{externalBranchSplitFnc[[1]]}.}

  \item{pamStage}{Only used for method "hybrid".  If TRUE, the second (PAM-like) stage will be performed.  }  

  \item{pamRespectsDendro}{Logical, only used for method "hybrid". If \code{TRUE}, the PAM stage will
respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters
that belong to the same branch that the objects or small clusters being assigned belong to.
}    

  \item{useMedoids}{Only used for method "hybrid" and only if \code{labelUnlabeled==TRUE}.
If TRUE, the second stage will be use object to medoid distance; if FALSE, it
  will use average object to cluster distance. The default (FALSE) is recommended. }

  \item{maxDistToLabel}{Deprecated, use \code{maxPamDist} instead. Only used for method "hybrid" and only if
\code{labelUnlabeled==TRUE}.
Maximum object distance to closest cluster that will result in the object
assigned to that cluster. }

  \item{maxPamDist}{Only used for method "hybrid" and only if \code{labelUnlabeled==TRUE}. 
Maximum object distance to closest cluster that will result in the object
assigned to that cluster. Defaults to \code{cutHeight}. }

  \item{respectSmallClusters}{Only used for method "hybrid" and only if \code{labelUnlabeled==TRUE}.
If TRUE, branches that failed to be clusters in stage 1 only because of
insufficient size will be assigned together in stage 2. If FALSE, all objects will be assigned
individually.}

  \item{verbose}{Controls the verbosity of the output. 0 will make the function completely quiet,
values up to 4 gradually increase verbosity.}

  \item{indent}{Controls indentation of printed messages (see \code{verbose} above). Each unit
adds two spaces before printed messages; useful when several functions' output is to be nested. }

}
\details{
  This is a wrapper for two related but different methods for cluster detection in hierarchical
clustering dendrograms.

In order to make the shape parameters \code{maxCoreScatter} and \code{minGap} more universal, their
values are interpreted relative to \code{cutHeight} and the 5th percetile of the merging heights (we
arbitrarily chose the 5th percetile rather than the minimum for reasons of stability). Thus, the absolute
maximum allowable core scatter is calculated as \code{maxCoreScatter * (cutHeight - refHeight) +
refHeight} and the absolute minimum allowable gap as \code{minGap * (cutHeight - refHeight)}, where
\code{refHeight} is the 5th percentile of the merging heights.
}
\value{
  A vector of numerical labels giving assignment of objects to modules. Unassigned objects are labeled
0, the largest module has label 1, next largest 2 etc.
}
\references{
  Langfelder P, Zhang B, Horvath S, 2007.
  \url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/BranchCutting}
}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{
  \code{\link{hclust}}, \code{\link{cutreeHybrid}}, \code{\link{cutreeDynamicTree}}.
}
\keyword{cluster}
