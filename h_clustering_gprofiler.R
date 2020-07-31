
h_clustering_gprofiler<-function(signpat,kmin=2,kmax=4){
  # calculate binomial distance between samples
  distance <- vegan::vegdist(t(signpat), method="binomial")
  
  # cluster tree
  cluster <- stats::hclust(distance, method="complete") # hclust from "stats"
  
  # cut cluster based on k
  groups <- list()
  cnt <- 1
  for(k in kmin:kmax){
    groups[[cnt]] <- stats::cutree(cluster, k=k)
    cnt <- cnt+1
  }
  names(groups) <- kmin:kmax
  
  return(groups)
}