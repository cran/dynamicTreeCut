#----------------------------------------------------------------------------------------------
#
# cutreeDynamic
#
#----------------------------------------------------------------------------------------------
# A wrapper function for cutreeHybrid and cutreeDynamicTree.

cutreeDynamic = function(dendro, cutHeight = NULL, minClusterSize = 20, 
                       method = "hybrid", distM = NULL, deepSplit = (ifelse(method=="hybrid", 1, FALSE)), 
                       maxCoreScatter = NULL, minGap = NULL,
                       maxAbsCoreScatter = NULL, minAbsGap = NULL, 
                       pamStage = TRUE, pamRespectsDendro = TRUE,
                       useMedoids = FALSE, maxDistToLabel = cutHeight,
                       respectSmallClusters = TRUE, 
                       verbose = 2, indent = 0)
{

  #if (!is.null(labelUnlabeled))
  #{
  #  pamStage = labelUnlabeled;
  #  warning("The argument 'labelUnlabeled' is deprecated. Please use 'pamStage' instead.");
  #}

  if (class(dendro)!="hclust") stop("Argument dendro must have class hclust.");
  methods = c("hybrid", "tree");
  met = charmatch(method, methods);
  if ( (met==1) && (is.null(distM)) )
  {
    warning('cutreeDynamic: method "hybrid" requires a valid dissimilarity matrix "distM". Defaulting to method "tree".');
    met = 2;
  }
  if (is.na(met))
  {
    stop(paste("Invalid method argument. Accepted values are (unique abbreviations of)", 
                paste(methods, collapse = ", ")));
  } else if (met==1)
  {
    # if (is.null(distM)) stop('distM must be given when using method "hybrid"');
    return(cutreeHybrid(dendro = dendro, distM = distM, minClusterSize = minClusterSize, 
                      cutHeight = cutHeight, deepSplit = deepSplit,
                      maxCoreScatter = maxCoreScatter, minGap = minGap,
                      maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                      pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                      useMedoids = useMedoids, 
                      maxDistToLabel = maxDistToLabel, 
                      respectSmallClusters = respectSmallClusters, 
                      verbose = verbose, indent = indent)$labels);
  } else
  {
    return(cutreeDynamicTree(dendro = dendro, maxTreeHeight = cutHeight, deepSplit = deepSplit,
                             minModuleSize = minClusterSize)); 
  }
}
    
