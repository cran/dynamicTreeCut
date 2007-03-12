# Tree cut

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

.NewBranch = function(IsBasic = TRUE, IsClosed = FALSE,
            Content = NULL, MergingHeights = NULL,
            RootHeight = 1, LastMerge = 0, Size = 0, MergedInto = 0, Singletons = NULL,
            SingletonHeights = NULL, 
            AttachHeight = NULL, FailSize = FALSE )
{
  list(IsBasic = IsBasic, IsClosed = IsClosed,
       Content = Content, MergingHeights = MergingHeights,
       RootHeight = RootHeight, LastMerge = LastMerge, Size = Size, MergedInto = MergedInto,
       Singletons = Singletons, SingletonHeights = SingletonHeights,
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

cutreeHybrid = function(dendro, distM, cutHeight = NULL, minClusterSize = 20, deepSplit = 1,
                       maxCoreScatter = NULL, minGap = NULL, 
                       maxAbsCoreScatter = NULL, minAbsGap = NULL, clusterTrim = 0,
                       labelUnlabeled = TRUE,
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

  if ( (dim(distM)[1] != nMerge+1) | (dim(distM)[2]!=nMerge+1) ) 
    stop("distM has incorrect dimensions.");

  diag(distM) = 0;

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

  defMCS = c(0.64, 0.73, 0.82, 0.91);
  defMG = (1-defMCS)*3/4;
  nSplitDefaults = length(defMCS);

  # Convert deep split to range 1..4
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
    
  # For each merge, record the cluster that it belongs to
  IndMergeToBranch = rep(0, times = nMerge)

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
             .NewBranch(Content = dendro$merge[merge,], MergingHeights = rep(dendro$height[merge], 2), 
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
      Branches[[clust]]$Content = c(Branches[[clust]]$Content, gene);
      if (Branches[[clust]]$IsBasic) Branches[[clust]]$Singletons = c(Branches[[clust]]$Singletons, -gene);
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
        Branches[[smaller]]$IsClosed = TRUE;
        Branches[[smaller]]$FailSize = SmallerFailSize;
        Branches[[smaller]]$MergedInto = larger; 
        Branches[[smaller]]$AttachHeight = dendro$height[merge];
        Branches[[larger]]$Content = c(Branches[[larger]]$Content, smaller);
        if (Branches[[larger]]$IsBasic) 
        {
           Branches[[larger]]$Singletons = 
                   c(Branches[[larger]]$Singletons, Branches[[smaller]]$Singletons);
           Branches[[larger]]$SingletonHeights = 
                   c(Branches[[larger]]$SingletonHeights, Branches[[smaller]]$SingletonHeights);
        }
        Branches[[larger]]$MergingHeights = c(Branches[[larger]]$MergingHeights, dendro$height[merge]);
        Branches[[larger]]$Size = Branches[[larger]]$Size + Branches[[smaller]]$Size;
        Branches[[larger]]$LastMerge = merge;
        IndMergeToBranch[merge] = larger;
        RootBranch = larger;
      } else
      {
        # Close both clusters and start a composite cluster.
        nBranches = nBranches + 1;
        Branches[[smaller]]$IsClosed = TRUE;
        Branches[[larger]]$IsClosed = TRUE;
        Branches[[smaller]]$AttachHeight = dendro$height[merge];
        Branches[[larger]]$AttachHeight = dendro$height[merge];
        # print(paste("  Starting a composite cluster with number", nBranches));
        Branches[[nBranches]] = .NewBranch(IsBasic = FALSE, Content = clusts, 
                        MergingHeights = rep(dendro$height[merge], 2), Size = sum(sizes),
                        LastMerge = merge, 
                        #Singletons = c(Branches[[larger]]$Singletons, Branches[[smaller]]$Singletons)
                        Singletons = NULL
                        );
        IndMergeToBranch[merge] = nBranches;
        RootBranch = nBranches;
      }
    }
  }

  if (verbose>2) printFlush(paste(spaces, "..Going through detected branches and marking clusters.."));
  
  nPoints = nMerge+1;

  IsBasic = rep(TRUE, times = nBranches);
  IsBranch = rep(FALSE, times = nBranches);
  SmallLabels = rep(0, times = nPoints);
  Trimmed = rep(0, times = nPoints);
  for (clust in 1:nBranches)
  {
    if (is.null(Branches[[clust]]$AttachHeight)) Branches[[clust]]$AttachHeight = cutHeight;
    IsBasic[clust] = Branches[[clust]]$IsBasic;
    if (Branches[[clust]]$IsBasic)
    {
       coresize = .CoreSize(length(Branches[[clust]]$Singletons), minClusterSize);
       Core = Branches[[clust]]$Singletons[c(1:coresize)];
       Branches[[clust]]$Core = Core;
       CoreScatter = mean(apply(distM[Core, Core], 2, sum)/(coresize-1)); 
    } else { CoreScatter = 0; }
    IsBranch[clust] = Branches[[clust]]$IsBasic & (Branches[[clust]]$Size >= minClusterSize) & 
                (CoreScatter < maxAbsCoreScatter) &
                (Branches[[clust]]$AttachHeight - CoreScatter > minAbsGap);
    if (Branches[[clust]]$FailSize) SmallLabels[Branches[[clust]]$Singletons] = clust;
  }
  if (!respectSmallClusters) SmallLabels = rep(0, times = nPoints);

  # Trim objects from clusters that are too close to the boundary.

  if (clusterTrim>0) for (clust in 1:nBranches) 
    if (Branches[[clust]]$IsBasic & (Branches[[clust]]$Size>2))
  {
    if (is.null(Branches[[clust]]$AttachHeight)) Branches[[clust]]$AttachHeight = cutHeight;
    #if (length(Branches[[clust]]$SingletonHeights)!=length(Branches[[clust]]$Singletons))
    #  stop("Internal error: length of SingletonHeights differs from length of Singletons. Sorry!");
    bottom = min(Branches[[clust]]$MergingHeights);
    top = Branches[[clust]]$AttachHeight;
    # First object in the cluster to be trimmed:
    FirstTrim = match( TRUE, (Branches[[clust]]$MergingHeights - bottom)/(top-bottom) > 1-clusterTrim);
    if (!is.na(FirstTrim))
    {
      FirstCont = Branches[[clust]]$Content[FirstTrim[1]];
      NSingls = Branches[[clust]]$Size;
      if (FirstCont<0)
      {
        FirstSingl = c(1:NSingls)[Branches[[clust]]$Singletons == -FirstCont];
        if (length(FirstSingl)==0) 
        { 
           print(paste("FirstCont:", FirstCont))
           print("Content:");
           print( Branches[[clust]]$Content);
           print(paste("FirstSingl:", FirstSingl))
           print("Singletons:");
           print( Branches[[clust]]$Singletons);
           stop(paste("Internal error: Trimming: First trimmed content",
                                      "points to an invalid singleton. Sorry!"));
        }
        
      } else
      {
        FirstSingl = 
             c(1:NSingls)[Branches[[clust]]$Singletons == Branches[[FirstCont]]$Singletons[1]];
        if (length(FirstSingl)==0) stop(paste("Internal error: Trimming: First trimmed content",
                                      "points to an invalid cluster singleton. Sorry!"));
      }
      if ((FirstSingl < NSingls/2) | (FirstSingl < 3))
      {
        #printFlush(paste("Trimming cluster ", clust));
        #printFlush("Merging heights:")
        #print(Branches[[clust]]$MergingHeights);
        #printFlush(paste("FirstTrim:", FirstTrim));
        #printFlush("Trim Condition:");
        #printFlush((Branches[[clust]]$MergingHeights - bottom)/(top-bottom) <= clusterTrim);
        warning(paste("GetClusters: keeping a low proportion of a cluster:", FirstSingl,
                      "out of", Branches[[clust]]$Size, "objects.\n", 
                      "Increasing the proportion of kept objects to preserve a branch."));
        if (verbose>3)
        {
          printFlush(paste("GetClusters: keeping a low proportion of a cluster:", FirstSingl,
                      "out of", Branches[[clust]]$Size, "objects."));
          printFlush(paste(spaces, "  "));
          printFlush(paste(spaces, "  ..SingletonHeights:", 
                            paste(signif(Branches[[clust]]$SingletonHeights,2), collapse=", ")));
          printFlush(paste(spaces, "  ..bottom =", signif(bottom,2), ", top =", signif(top,2)));
        }
        # Make sure we keep a certain minimum of singletons
        FirstSingl = max(FirstSingl, as.integer(NSingls/3)+1, 3);
      }
      if (FirstSingl<=NSingls) 
      {
        Trimmed[Branches[[clust]]$Singletons[c(FirstSingl:NSingls)]] = clust;
        Branches[[clust]]$Singletons = Branches[[clust]]$Singletons[-c(FirstSingl:NSingls)];
        Branches[[clust]]$SingletonHeights = Branches[[clust]]$SingletonHeights[-c(FirstSingl:NSingls)];
        Branches[[clust]]$Size = length(Branches[[clust]]$Singletons);
      }
    }
  }

  # Here's where the original AssignLabel starts
  
  if (verbose>2) printFlush(paste(spaces, "..Assigning stage 1 labels.."));

  Colors = rep(0, times = nPoints);
  IsCore = rep(0, times = nPoints);
  BranchBranches = c(1:nBranches)[IsBranch];
  color = 0;
  for (clust in BranchBranches)
  {
    color = color+1;
    Colors[Branches[[clust]]$Singletons] = color;
    SmallLabels[Branches[[clust]]$Singletons] = 0;
    coresize = .CoreSize(length(Branches[[clust]]$Singletons), minClusterSize);
    Core = Branches[[clust]]$Singletons[c(1:coresize)];
    IsCore[Core] = color;
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
  if (labelUnlabeled & UnlabeledExist & nProperLabels>0)
  {
     if (verbose>1) printFlush(paste(spaces, "..Assigning stage 2 labels.."));
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
            # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
            DistInCluster = distM[InCluster, InCluster];
            if (length(InCluster)>1)
            {
              DistSums = apply(DistInCluster, 2, sum);
              smed = InCluster[which.min(DistSums)];
              DistToMeds = distM[Medoids, smed];
              Nearest = which.min(DistToMeds);
              DistToNearest = DistToMeds[Nearest];
              if ( (DistToNearest < ClusterRadii[Nearest]) | (DistToNearest <  maxDistToLabel) )
              {
                Colors[InCluster] = Nearest;
              } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later 
            }
          }
        }
        # Assign leftover unlabeled objects to clusters with nearest medoids
        Unlabeled = c(1:nPoints)[Colors==0];
        UnassdToMedoidDist = distM[Medoids, Unlabeled];
        if (nProperLabels>1)
        {
           NearestMedoids = apply(UnassdToMedoidDist, 2, which.min);
           NearestCenterDist = apply(UnassdToMedoidDist, 2, min);
        } else {
           NearestMedoids = rep(1, times = nUnlabeled);
           NearestCenterDist = UnassdToMedoidDist;
        }
        Colors[Unlabeled] = ifelse((NearestCenterDist < ClusterRadii[NearestMedoids]) | 
                                    (NearestCenterDist < maxDistToLabel) ,
                                    NearestMedoids, 0);
        UnlabeledExist = (sum(Colors==0)>0);
        if (verbose>1)
          printFlush(paste(spaces, "   ...assigned", 
                          sum((NearestCenterDist < ClusterRadii[NearestMedoids]) | 
                                    (NearestCenterDist < maxDistToLabel)), 
                            "of", nUnlabeled, "previously unassigned points.")); 
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
          if (nSmallClusters>0) for (sclust in SmallLabLevs[SmallLabLevs!=0])
          {
            InCluster = c(1:nPoints)[SmallLabels==sclust];
            # printFlush(paste("SmallCluster", sclust, "has", length(InCluster), "elements."));
            if (length(InCluster)>1)
            {
              DistSClustClust = distM[InCluster, Labeled];
              MeanDist = apply(DistSClustClust, 2, mean);
              MeanMeanDist = tapply(MeanDist, LabelFac, mean);
              Nearest = which.min(MeanMeanDist);
              NearestDist = MeanMeanDist[Nearest];
              if ( ((NearestDist < ClusterDiam[Nearest]) | (NearestDist <  maxDistToLabel)) )
              {
                Colors[InCluster] = Nearest;
              } else Colors[InCluster] = -1;  # This prevents individual points from being assigned later
            } 
          }
        }
        # Assign leftover unlabeled objects to clusters with nearest medoids
        Unlabeled = c(1:nPoints)[Colors==0];
        if (length(Unlabeled)>0)
        {
          if (length(Unlabeled)>1)
          {
            UnassdToClustDist = apply(distM[Labeled, Unlabeled], 2, tapply, LabelFac, mean);
            if (nProperLabels>1)
            {
               NearestClusters = apply(UnassdToClustDist, 2, which.min);
               NearestClusterDist = apply(UnassdToClustDist, 2, min);
            } else
            {
               NearestClusters = rep(1, length(Unlabeled));
               NearestClusterDist = UnassdToClustDist;
            }
          } else {
            UnassdToClustDist = tapply(distM[Labeled, Unlabeled], LabelFac, mean);
            NearestClusters = which.min(UnassdToClustDist);
            NearestClusterDist = min(UnassdToClustDist);
          }
          Colors[Unlabeled] = ifelse((NearestClusterDist < ClusterDiam[NearestClusters]) | 
                                      (NearestClusterDist < maxDistToLabel),
                                       NearestClusters, 0);
          if (verbose>1)
            printFlush(paste(spaces, "   ...assigned", 
                            sum((NearestClusterDist < ClusterDiam[NearestClusters]) | 
                                      (NearestClusterDist < maxDistToLabel)), 
                              "of", nUnlabeled, "previously unassigned points.")); 
        }
     }
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
       trimmed = as.numeric(factor(Trimmed))-1,
       branches  = list(nBranches = nBranches, Branches = Branches, 
                        IndMergeToBranch = IndMergeToBranch,
                        RootBranch = RootBranch, IsBasic = IsBasic, IsBranch = IsBranch, 
                        nPoints = nMerge+1));
} 


#----------------------------------------------------------------------------------------------
#
# cutreeDynamic
#
#----------------------------------------------------------------------------------------------
# A wrapper function for cutreeHybrid and cutreeDynamicTree.

cutreeDynamic = function(dendro, cutHeight = NULL, minClusterSize = 20, 
                       method = "hybrid", distM = NULL, deepSplit = (ifelse(method=="hybrid", 1, FALSE)), 
                       maxCoreScatter = NULL, minGap = NULL,
                       maxAbsCoreScatter = NULL, minAbsGap = NULL, clusterTrim = 0,  
                       labelUnlabeled = TRUE,
                       useMedoids = FALSE, maxDistToLabel = cutHeight,
                       respectSmallClusters = TRUE, 
                       verbose = 2, indent = 0)
{

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
    return(cutreeHybrid(dendro = dendro, distM = as.matrix(distM), minClusterSize = minClusterSize, 
                      cutHeight = cutHeight, deepSplit = deepSplit,
                      maxCoreScatter = maxCoreScatter, minGap = minGap,
                      maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                      labelUnlabeled = labelUnlabeled, useMedoids = useMedoids, 
                      maxDistToLabel = maxDistToLabel, clusterTrim = clusterTrim, 
                      respectSmallClusters = respectSmallClusters, 
                      verbose = verbose, indent = indent)$labels);
  } else
  {
    return(cutreeDynamicTree(dendro = dendro, maxTreeHeight = cutHeight, deepSplit = deepSplit,
                             minModuleSize = minClusterSize)); 
  }
}
    
