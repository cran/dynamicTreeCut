
# Keep the name dynamicTreeCut; argument expr can be NULL; if NULL, just revert to the standard tree cut
# that does not look at the second PC.

# 1.20:

# Another improvement that would be beneficial is to optionally restrict PAM assignment to only the
# branch on which the gene sits.

# ClusterTrim is removed for now; it would have been a bit more complicated to make it work with the
# enhanced PAM stage. May be re-introduced later if anyone shows any interest in this.

# 1.11-2
#    . Bug fix in interpretation of deepSplit fixed

# 1.11-1
#    . If no merge lies below the cut height, simply exit with all labels=0 instead of throwing an
#    error.
#    . If cutHeight is above the highest merge, set it to the highest merge.

# 1.11
#    . change the default cutHeight to 99% of dendogram height range
 
# 1.10-02:
#    . if distM isn't given, method defaults to "tree".

# 1.10:
#    . Bug fix: since distM us necessary in the first stage as well, the function now complains if
#    distM is not given

# 1.09:
#    . Bug fix: number of Unassigned = 1 is now handled correctly

# 1.08: 
#    . Changing the meaning of minGap, maxCoreScatter: both now relative (i.e., fractions). minGap will
#    be the minimum cluster gap expressed as a fraction of the range between a certain quantile of the
#    merging heights and the cutHeight; maxCoreScatter will be interpreted in the same way. Adding two
#    more parameters, namely minAbsGap, maxAbsCoreScatter that can be used to supply hitherto used
#    absolute values. If they are given, they override minGap and maxCoreScatter, respectively. 
#    . Change of names: clusterMinSize -> minClusterSize
#    . Default value of minClusterSize now 20
#    . Fixed stage 2 labeling for non-medoids: ClusterDiam{eter} now given is the maximum of average
#      distances of points to the rest of the cluster.
#    . Default for maxDistToLabel is now cutHeight even though it may mean that some objects above the
#      cutHeight may be labeled. There doesn't seem to be a clean and simple way to only label objects
#      whose merging distance to (the continuation of) the cluster is below cutHeight
#    . If cutHeight not given, will be set to max(dendro$height) for the DynamicTree as well.

# 1.07: Changing parameter names to be more intuitive and in accord with generally used terminology
#       Also changing internal function names by prepending a . to them
#    . All extra functions (EvaluateCLusters etc) removed.
#    . Rename the function from GetClusters to cutreeHybrid

# 1.06: The tree cut is exactly the same as 1.05, but the color handling is relegated to NetworkFunctions.
#       The ColorsFromLabels in NetworkFunctions use a slightly different input format; 
#           in particular, 0 is considered grey. This means, among other things, that 
# 1.05:
#   . The PAM stage is changed: instead of calculating distances to medoids, will calculate average
#     distances to clusters. This is intuitively less desirable than the medoids, but simulations
#     seem to indicate that large clusters are too greedy. The average linkage may help that a bit.
#     It is not quite clear that cluster trimming makes a lot of sense in this scenario, but I'll
#     keep it in (assigning elements based on average distance to a trimmed cluster is not quite the
#     same as an untrimmed cluster, even for the elements that were trimmed). 
#   . Improve trimming: instead of the lowest joining heights, keep all singletons up to the first joined
#     object (branch or singleton) whose merging height is above the threshold; the rest (content of all
#     branches merged higher than threshold) is trimmed.

# 1.04:
#   . Implement Bin's idea that small branches unassigned in Stage 1 should not be broken
#     up by Stage 2. Implemented as follows: Only keep those branches unbroken whose only reason for not
#     being a cluster is that they don't have enough objects.
#     While going through the merge tree: mark branches that are (1) merged into composite clusters, (2)
#     not clusters themselves because of failing the minimum size requirement only.
#   . Change boudaries of clusters. Instead of regarding everything up to the merge automatically part
#     the respective cluster, only elements whose joining heights are less than a cutoff given for
#     the cluster are considered part of the cluster automatically; the rest is assigned in pam-like
#     manner. 

# 1.03: 
#   . Changing the definition of the core from the most connected points to the first points added
#     to the cluster. Note that the order depends on how branches are merged, so that better be
#     correct as well. This makes the core stable against adding outliers to the cluster.

# 1.02.01: fixing a bug in the main function that was referencing nonexistent heights member of
#          dendro.
#    . Fixing a correctness issue: the core average distance is the average distance between points
#      within the core, not the average distance of points within the core to all points in the
#      cluster. 

# 1.02: 
#    . Changing core size: instead of minClusterSize + sqrt(Size - minClusterSize) it is now
#      CoreSize  + sqrt(Size - CoreSize), with CoreSize = as.integer(minClusterSize) + 1.


# 1.01.03: Fixing memory usage and a bug in which singletons were added twice (and more times).
#  In this version, to simplify things, Singletons are only kept for basic clusters. For composite
#  clusters they are NULL; CoreScatter should never be called for composite clusters as that can get
#  hugely time and memory intensive.
# Another bug concerning couting unassigned points fixed.

# -02: adding distance matrix to the parameters of GetClusters.
# -03: Merging GetClusters and AssignLabel together; changes in variable names to make code more
# readable.

.NewBranch = function(IsBasic = TRUE, IsTopBasic = IsBasic,
            Singletons = NULL,
            Clusters = NULL,
            BasicClusters = NULL,
            MergingHeights = NULL,
            RootHeight = 1, LastMerge = 0, Size = 0, MergedInto = 0,
            SingletonHeights = NULL, 
            AttachHeight = NULL, FailSize = FALSE )
{
  list(IsBasic = IsBasic, 
       IsTopBasic = IsTopBasic,
       Singletons = Singletons,
       Clusters = Clusters,
       BasicClusters = BasicClusters,
       MergingHeights = MergingHeights,
       RootHeight = RootHeight, LastMerge = LastMerge, Size = Size, MergedInto = MergedInto,
       SingletonHeights = SingletonHeights,
       AttachHeight = AttachHeight, FailSize = FailSize);
}

# The following are supporting function for GetClusters. 

.CoreSize = function(BranchSize, minClusterSize)
{
  BaseCoreSize = minClusterSize/2 + 1;
  if (BaseCoreSize < BranchSize)
  {
    CoreSize = as.integer(BaseCoreSize + sqrt(BranchSize - BaseCoreSize));
  } else CoreSize = BranchSize;
  CoreSize;
}
  
# This assumes the diagonal of the distance matrix
# is zero, BranchDist is a square matrix whose dimension is at least 2.

.CoreScatter = function(BranchDist, minClusterSize)
{
  nPoints = dim(BranchDist)[1];
  PointAverageDistances = apply(BranchDist, 2, sum) / (nPoints-1);
  CoreSize = minClusterSize/2 + 1;
  if (CoreSize < nPoints)
  {
    EffCoreSize = as.integer(CoreSize + sqrt(nPoints - CoreSize));
    ord = order(PointAverageDistances);
    Core = ord[c(1:EffCoreSize)];
  } else {
    Core = c(1:nPoints);
    EffCoreSize = nPoints;
  }
  CoreAverageDistances = apply(BranchDist[Core, Core], 2, sum) / (EffCoreSize-1);
  mean(CoreAverageDistances);
}

#-------------------------------------------------------------------------------------------
#
# cutreeHybrid
#
#-------------------------------------------------------------------------------------------
# Traverses a given clustering tree and detects branches whose size is at least minClusterSize, average
# singleton joining height is at most maxCoreScatter and split (attaching height minus average
# height) is at least minGap. If cutHeight is set, all clusters are cut at that height.

# clusterTrim is the fraction of the cluster gap that will be trimmed away.
# Objects whose joining height is above that will be (re-)assigned based on distance to medoids. If
# clusterTrim<=0, all assigments of stage 1 will be respected. 

cutreeHybrid = function(dendro, distM, 
                       cutHeight = NULL, minClusterSize = 20, deepSplit = 1,
                       maxCoreScatter = NULL, minGap = NULL, 
                       maxAbsCoreScatter = NULL, minAbsGap = NULL,
                       pamStage = TRUE, pamRespectsDendro = TRUE,
                       useMedoids = FALSE, maxDistToLabel = cutHeight, 
                       respectSmallClusters = TRUE, 
                       verbose = 2, indent = 0)
{

  spaces = indentSpaces(indent);

  MxBranches = length(dendro$height)
  Branches = vector(mode="list", length = MxBranches);

  nBranches = 0;

  if (verbose>0) printFlush(paste(spaces, "Detecting clusters..."));

  # No. of merges in the tree
  nMerge = length(dendro$height);
  if (nMerge < 1)
    stop("The given dendrogram is suspicious: number of merges is zero.");

  if (is.null(distM)) stop("distM must be non-NULL")

  if (is.null(dim(distM)))
    stop("distM must be a matrix.");

  if ( (dim(distM)[1] != nMerge+1) | (dim(distM)[2]!=nMerge+1) ) 
    stop("distM has incorrect dimensions.");

  if (sum(abs(diag(distM)))>0) diag(distM) = 0;

  refQuantile = 0.05;
  refMerge = as.integer(nMerge * refQuantile + 0.5);
  if (refMerge < 1) refMerge = 1;
  refHeight = dendro$height[refMerge]; 
  
  if (is.null(cutHeight))
  {
    cutHeight = 0.99 * (max(dendro$height) - refHeight) + refHeight;
    if (verbose>0) 
      printFlush(paste(spaces, 
            "..cutHeight not given, setting it to", signif(cutHeight,3), 
            " ===>  99% of the (truncated) height range in dendro."));
  } else {
    if (cutHeight > max(dendro$height)) cutHeight = max(dendro$height);
  }

  nMergeBelowCut = sum(dendro$height <= cutHeight);
  if (nMergeBelowCut < minClusterSize) 
  {
    if (verbose>0) printFlush(paste(spaces, "cutHeight set too low: no merges below the cut."));
    return(list(labels = rep(0, times = nMerge+1)))
  }

  # Default values for maxCoreScatter and minGap:

  defMCS = c(0.64, 0.73, 0.82, 0.91, 0.95);
  defMG = (1-defMCS)*3/4;
  nSplitDefaults = length(defMCS);

  # Convert deep split to range 1..5
  if (is.logical(deepSplit)) deepSplit = as.integer(deepSplit)*(nSplitDefaults - 2);
  deepSplit = as.integer(deepSplit + 1);
  if ((deepSplit<1) | (deepSplit>nSplitDefaults))
    stop(paste("Parameter deepSplit (value", deepSplit, 
               ") out of range: allowable range is 0 through", nSplitDefaults-1));

  # If not set, set the cluster gap and core scatter according to deepSplit.
  if (is.null(maxCoreScatter)) maxCoreScatter = defMCS[deepSplit];
  if (is.null(minGap)) minGap = defMG[deepSplit];

  # If maxDistToLabel is not set, set it to cutHeight

  if (is.null(maxDistToLabel)) maxDistToLabel = cutHeight;

  # Convert (relative) minGap and maxCoreScatter to corresponding absolute quantities if the latter were
  # not given.

  if (is.null(maxAbsCoreScatter))
     maxAbsCoreScatter = refHeight + maxCoreScatter * (cutHeight - refHeight);
  if (is.null(minAbsGap))
     minAbsGap = minGap * (cutHeight - refHeight);
    
  nPoints = nMerge+1;

  # For each merge, record the cluster that it belongs to
  IndMergeToBranch = rep(0, times = nMerge)

  # For each object that joins a composite branch, record the number of the branch
  onBranch = rep(0, nPoints);

  # The root
  RootBranch = 0;
  if (verbose>2) printFlush(paste(spaces, "..Going through the merge tree"));

  for (merge in 1:nMerge) if (dendro$height[merge]<=cutHeight)
  {
    # are both merged objects sigletons?
    if (dendro$merge[merge,1]<0 & dendro$merge[merge,2]<0)
    {
      # Yes; start a new cluster.
      nBranches = nBranches + 1;
      Branches[[nBranches]] = 
             .NewBranch(MergingHeights = rep(dendro$height[merge], 2), 
                        LastMerge = merge, Size = 2, Singletons = -dendro$merge[merge,],
                        SingletonHeights = rep(dendro$height[merge], 2));
      IndMergeToBranch[merge] = nBranches;
      RootBranch = nBranches;
    } else if (dendro$merge[merge,1] * dendro$merge[merge,2] <0)
    {
      # merge the sigleton into the cluster
      clust = IndMergeToBranch[max(dendro$merge[merge,])];
      if (clust==0) stop("Internal error: a previous merge has no associated cluster. Sorry!");
      gene = min(dendro$merge[merge,]);
      if (Branches[[clust]]$IsBasic) 
      {
        Branches[[clust]]$Singletons = c(Branches[[clust]]$Singletons, -gene);
      } else {
        onBranch[-gene] = clust;
      }
      Branches[[clust]]$MergingHeights = c(Branches[[clust]]$MergingHeights, dendro$height[merge]);
      Branches[[clust]]$SingletonHeights = c(Branches[[clust]]$SingletonHeights, dendro$height[merge]);
      Branches[[clust]]$LastMerge = merge;
      Branches[[clust]]$Size = Branches[[clust]]$Size + 1;
      IndMergeToBranch[merge] = clust;
      RootBranch = clust;
    } else
    {
      # attempt to merge two clusters:
      clusts = IndMergeToBranch[dendro$merge[merge,]];
      sizes = c(Branches[[clusts[1]]]$Size, Branches[[clusts[2]]]$Size);
      rnk = rank(sizes, ties.method = "first");
      smaller = clusts[rnk[1]]; larger = clusts[rnk[2]];
      if (Branches[[smaller]]$IsBasic)
      {
         coresize = .CoreSize(length(Branches[[smaller]]$Singletons), minClusterSize);
         Core = Branches[[smaller]]$Singletons[c(1:coresize)];
         SmAveDist = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
      } else { SmAveDist = 0; }
       
      if (Branches[[larger]]$IsBasic)
      {
         coresize = .CoreSize(length(Branches[[larger]]$Singletons), minClusterSize);
         Core = Branches[[larger]]$Singletons[c(1:coresize)];
         LgAveDist = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
      } else { LgAveDist = 0; }

      # Is the smaller cluster small or shallow enough to be merged?
      SmallerScores = c(Branches[[smaller]]$IsBasic, 
                        (Branches[[smaller]]$Size < minClusterSize),
                        (SmAveDist > maxAbsCoreScatter), 
                        (dendro$height[merge] - SmAveDist < minAbsGap));
      
      if ( SmallerScores[1] * sum(SmallerScores[c(2:4)]) > 0 )
      {
        DoMerge = TRUE;
        SmallerFailSize = !(SmallerScores[3] | SmallerScores[4]);  # Smaller fails only due to size
      } else 
      {
        LargerScores = c(Branches[[larger]]$IsBasic, 
                          (Branches[[larger]]$Size < minClusterSize),
                          (LgAveDist > maxAbsCoreScatter), 
                          (dendro$height[merge] - LgAveDist < minAbsGap));
        if ( LargerScores[1] * sum(LargerScores[c(2:4)]) > 0 )
        { # Actually: the larger one is the one to be merged
          DoMerge = TRUE;
          SmallerFailSize = !(LargerScores[3] | LargerScores[4]);  # cluster fails only due to size
          x = smaller; smaller = larger; larger = x;
        } else {
          DoMerge = FALSE; # None of the two satisfies merging criteria
        }
      }
      if (DoMerge)
      {
        # merge the smaller into the larger cluster and close it.
        Branches[[smaller]]$FailSize = SmallerFailSize;
        Branches[[smaller]]$MergedInto = larger; 
        Branches[[smaller]]$AttachHeight = dendro$height[merge];
        if (Branches[[larger]]$IsBasic) 
        {
           Branches[[larger]]$Singletons = 
                   c(Branches[[larger]]$Singletons, Branches[[smaller]]$Singletons);
           Branches[[larger]]$SingletonHeights = 
                   c(Branches[[larger]]$SingletonHeights, Branches[[smaller]]$SingletonHeights);
           Branches[[smaller]]$IsTopBasic = FALSE;
        } else {
          if (!Branches[[smaller]]$IsBasic)
            stop("Internal error: merging two composite clusters. Sorry!");
          onBranch[Branches[[smaller]]$Singletons] = larger
        }
        Branches[[larger]]$MergingHeights = c(Branches[[larger]]$MergingHeights, dendro$height[merge]);
        Branches[[larger]]$Size = Branches[[larger]]$Size + Branches[[smaller]]$Size;
        Branches[[larger]]$LastMerge = merge;
        IndMergeToBranch[merge] = larger;
        RootBranch = larger;
      } else
      {
        # start a composite cluster.
        nBranches = nBranches + 1;
        Branches[[smaller]]$AttachHeight = dendro$height[merge];
        Branches[[larger]]$AttachHeight = dendro$height[merge];
        Branches[[smaller]]$MergedInto = nBranches;
        Branches[[larger]]$MergedInto = nBranches;
        if (Branches[[smaller]]$IsBasic)
        {
          contdBasicClusters = smaller; # contained basic clusters
        } else {
          contdBasicClusters = Branches[[smaller]]$BasicClusters;
        }
        if (Branches[[larger]]$IsBasic)
        {
          contdBasicClusters = c(contdBasicClusters, larger);
        } else {
          contdBasicClusters = c(contdBasicClusters, Branches[[larger]]$BasicClusters);
        }
        # print(paste("  Starting a composite cluster with number", nBranches));
        Branches[[nBranches]] = .NewBranch(IsBasic = FALSE, Clusters = clusts, 
                        BasicClusters = contdBasicClusters,
                        MergingHeights = rep(dendro$height[merge], 2), Size = sum(sizes),
                        LastMerge = merge, 
                        Singletons = NULL);
        IndMergeToBranch[merge] = nBranches;
        RootBranch = nBranches;
      }
    }
  }

  if (verbose>2) printFlush(paste(spaces, "..Going through detected branches and marking clusters.."));
  
  IsBasic = rep(TRUE, times = nBranches);
  IsBranch = rep(FALSE, times = nBranches);
  SmallLabels = rep(0, times = nPoints);
  for (clust in 1:nBranches)
  {
    if (is.null(Branches[[clust]]$AttachHeight)) Branches[[clust]]$AttachHeight = cutHeight;
    IsBasic[clust] = Branches[[clust]]$IsBasic;
    if (Branches[[clust]]$IsTopBasic)
    {
       coresize = .CoreSize(length(Branches[[clust]]$Singletons), minClusterSize);
       Core = Branches[[clust]]$Singletons[c(1:coresize)];
       Branches[[clust]]$Core = Core;
       CoreScatter = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
       IsBranch[clust] = Branches[[clust]]$IsTopBasic & (Branches[[clust]]$Size >= minClusterSize) & 
                  (CoreScatter < maxAbsCoreScatter) &
                  (Branches[[clust]]$AttachHeight - CoreScatter > minAbsGap);
    } else { CoreScatter = 0; }
    if (Branches[[clust]]$FailSize) SmallLabels[Branches[[clust]]$Singletons] = clust;
  }
  if (!respectSmallClusters) SmallLabels = rep(0, times = nPoints);

  if (verbose>2) printFlush(paste(spaces, "..Assigning Tree Cut stage labels.."));

  Colors = rep(0, times = nPoints);
  IsCore = rep(0, times = nPoints);
  BranchBranches = c(1:nBranches)[IsBranch];
  branchLabels = rep(0, length(Branches));
  color = 0;
  for (clust in BranchBranches)
  {
    color = color+1;
    Colors[Branches[[clust]]$Singletons] = color;
    SmallLabels[Branches[[clust]]$Singletons] = 0;
    coresize = .CoreSize(length(Branches[[clust]]$Singletons), minClusterSize);
    Core = Branches[[clust]]$Singletons[c(1:coresize)];
    IsCore[Core] = color;
    branchLabels[clust] = color;
  } 

  Labeled = c(1:nPoints)[Colors!=0];
  Unlabeled = c(1:nPoints)[Colors==0];
  nUnlabeled = length(Unlabeled);
  UnlabeledExist = (nUnlabeled>0);
  if (length(Labeled)>0)
  {
    LabelFac = factor(Colors[Labeled]);
    nProperLabels = nlevels(LabelFac);
  } else
  {
    nProperLabels = 0;
  }
  if (pamStage & UnlabeledExist & nProperLabels>0)
  {
     if (verbose>2) printFlush(paste(spaces, "..Assigning PAM stage labels.."));
     nPAMed = 0;
     # Assign some of the grey genes to the nearest module. Define nearest as the distance to the medoid,
     # that is the point in the cluster that has the lowest average distance to all other points in the
     # cluster. First get the medoids.
     if (useMedoids)
     {
        Medoids = rep(0, times = nProperLabels);
        ClusterRadii = rep(0, times = nProperLabels);	
        for (cluster in 1:nProperLabels)
        {
          InCluster = c(1:nPoints)[Colors==cluster];
          DistInCluster = distM[InCluster, InCluster];
          DistSums = apply(DistInCluster, 2, sum);
          Medoids[cluster] = InCluster[which.min(DistSums)];
          ClusterRadii[cluster] = max(DistInCluster[, which.min(DistSums)])
        }
        # If small clusters are to be respected, assign those first based on medoid-medoid distances.
        if (respectSmallClusters)
        {
          FSmallLabels = factor(SmallLabels);
          SmallLabLevs = as.numeric(levels(FSmallLabels));
          nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1]==0);
          if (nSmallClusters>0) for (sclust in SmallLabLevs[SmallLabLevs!=0])
          {
            InCluster = c(1:nPoints)[SmallLabels==sclust];
            if (pamRespectsDendro)
            {
              onBr = unique(onBranch[InCluster]);
              if (length(onBr>1))
                stop(paste("Internal error: objects in a small cluster are marked to belong",
                           "\nto several large branches:", paste(onBr, collapse = ", ")));
              if (onBr > 0)
              {
                 basicOnBranch = Branches[[onBr]]$BasicClusters;
                 labelsOnBranch = branchLabels[basicOnBranch]
              } else {
                 labelsOnBranch = NULL;    
              }
            } else {
                labelsOnBranch = c(1:nProperLabels)
            }
            # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
            DistInCluster = distM[InCluster, InCluster, drop = FALSE];
            if (length(labelsOnBranch) > 0)
            {
               if (length(InCluster)>1)
               {
                 DistSums = apply(DistInCluster, 2, sum);
                 smed = InCluster[which.min(DistSums)];
                 DistToMeds = distM[Medoids[labelsOnBranch], smed];
                 closest = which.min(DistToMeds);
                 DistToClosest = DistToMeds[closest];
                 closestLabel = labelsOnBranch[closest]
                 if ( (DistToClosest < ClusterRadii[closestLabel]) | (DistToClosest <  maxDistToLabel) )
                 {
                   Colors[InCluster] = closestLabel;
                   nPAMed = nPAMed + length(InCluster);
                 } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later 
               }
            } else
              Colors[InCluster] = -1;
          }
        }
        # Assign leftover unlabeled objects to clusters with nearest medoids
        Unlabeled = c(1:nPoints)[Colors==0];
        if (length(Unlabeled>0)) for (obj in Unlabeled)
        {
          if (pamRespectsDendro)
          {
            onBr = onBranch[obj]
            if (onBr > 0)
            {
              basicOnBranch = Branches[[onBr]]$BasicClusters;
              labelsOnBranch = branchLabels[basicOnBranch]
            } else {
              labelsOnBranch = NULL;
            }
          } else {
              labelsOnBranch = c(1:nProperLabels)
          }
          if (!is.null(labelsOnBranch))
          {
             UnassdToMedoidDist = distM[Medoids[labelsOnBranch], obj];
             nearest= which.min(UnassdToMedoidDist)
             NearestCenterDist = UnassdToMedoidDist[nearest];
             nearestMed = labelsOnBranch[nearest]
             if ( (NearestCenterDist < ClusterRadii[nearestMed]) |
                  (NearestCenterDist < maxDistToLabel))
             {
               Colors[obj] = nearestMed;
               nPAMed = nPAMed + 1;
             }
          }
        }
        UnlabeledExist = (sum(Colors==0)>0);
     } else # Instead of medoids, use average distances
     {
        ClusterDiam = rep(0, times = nProperLabels);	
        for (cluster in 1:nProperLabels)
        {
          InCluster = c(1:nPoints)[Colors==cluster];
          nInCluster = length(InCluster)
          DistInCluster = distM[InCluster, InCluster];
          if (nInCluster>1) {
             AveDistInClust = apply(DistInCluster, 2, sum)/(nInCluster-1);
             ClusterDiam[cluster] = max(AveDistInClust);
          } else {
             ClusterDiam[cluster] = 0;
          }
        }
        # If small clusters are respected, assign them first based on average cluster-cluster distances.
        if (respectSmallClusters)
        {
          FSmallLabels = factor(SmallLabels);
          SmallLabLevs = as.numeric(levels(FSmallLabels));
          nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1]==0);
          ColorsX = Colors;
          if (nSmallClusters>0) for (sclust in SmallLabLevs[SmallLabLevs!=0])
          {
            InCluster = c(1:nPoints)[SmallLabels==sclust];
            if (pamRespectsDendro)
            {
              onBr = unique(onBranch[InCluster]);
              if (length(onBr)>1)
                stop(paste("Internal error: objects in a small cluster are marked to belong",
                           "\nto several large branches:", paste(onBr, collapse = ", ")));
              if (onBr > 0)
              {
                 basicOnBranch = Branches[[onBr]]$BasicClusters;
                 labelsOnBranch = branchLabels[basicOnBranch]
              } else {
                 labelsOnBranch = NULL;
              }
            } else {
                labelsOnBranch = c(1:nProperLabels)
            }
            # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
            if (!is.null(labelsOnBranch))
            {
              useObjects = is.finite(match(ColorsX, labelsOnBranch))
              DistSClustClust = distM[InCluster, useObjects, drop = FALSE];
              MeanDist = apply(DistSClustClust, 2, mean);
              useColorsFac = factor(ColorsX[useObjects])
              MeanMeanDist = tapply(MeanDist, useColorsFac, mean);
              nearest = which.min(MeanMeanDist);
              NearestDist = MeanMeanDist[nearest];
              nearestLabel = as.numeric(levels(useColorsFac)[nearest])
              if ( ((NearestDist < ClusterDiam[nearestLabel]) | (NearestDist <  maxDistToLabel)) )
              {
                Colors[InCluster] = nearestLabel;
                nPAMed = nPAMed + length(InCluster);
              } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later
            } 
          }
        }
        # Assign leftover unlabeled objects to clusters with nearest medoids
        Unlabeled = c(1:nPoints)[Colors==0];
        #ColorsX = Colors;
        if (length(Unlabeled)>0) for (obj in Unlabeled)
        {
          if (pamRespectsDendro)
          {
            onBr = onBranch[obj]
            if (onBr > 0)
            {
              basicOnBranch = Branches[[onBr]]$BasicClusters;
              labelsOnBranch = branchLabels[basicOnBranch]
            } else {
              labelsOnBranch = NULL;
            }
          } else {
            labelsOnBranch = c(1:nProperLabels)
          }
          if (!is.null(labelsOnBranch))
          {
            useObjects = is.finite(match(ColorsX, labelsOnBranch))
            useColorsFac = factor(ColorsX[useObjects])
            UnassdToClustDist = tapply(distM[useObjects, obj], useColorsFac, mean);
            nearest = which.min(UnassdToClustDist);
            NearestClusterDist = UnassdToClustDist[nearest];
            nearestLabel = as.numeric(levels(useColorsFac)[nearest])
            if ((NearestClusterDist < ClusterDiam[nearestLabel]) |
                (NearestClusterDist < maxDistToLabel) )
            {
              Colors[obj] = nearestLabel
              nPAMed = nPAMed + 1;
            }
          }
        }
     }
     if (verbose>2) printFlush(paste(spaces, "....assigned", nPAMed, "objects to existing clusters."));
  }

  # Relabel labels such that 1 corresponds to the largest cluster etc.
  Colors[Colors<0] = 0;
  UnlabeledExist = (sum(Colors==0)>0);
  NumLabs = as.numeric(as.factor(Colors));
  Sizes = table(NumLabs);
  if (UnlabeledExist)
  {
    if (length(Sizes)>1)
    {
        SizeRank = c(1, rank(-Sizes[2:length(Sizes)], ties.method="first")+1);
    } else {
        SizeRank = 1;
    }
    OrdNumLabs = SizeRank[NumLabs];
  } else {
    SizeRank = rank(-Sizes[1:length(Sizes)], ties.method="first");
    OrdNumLabs = SizeRank[NumLabs];
  }
  OrdIsCore = OrdNumLabs-UnlabeledExist;
  OrdIsCore[IsCore==0] = 0; 

  list(labels = OrdNumLabs-UnlabeledExist,
       cores = OrdIsCore,
       smallLabels = SmallLabels,
       onBranch = onBranch,
       branches  = list(nBranches = nBranches, Branches = Branches, 
                        IndMergeToBranch = IndMergeToBranch,
                        RootBranch = RootBranch, IsBasic = IsBasic, IsBranch = IsBranch, 
                        nPoints = nMerge+1));
} 


