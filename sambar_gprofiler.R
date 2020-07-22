sambar_gprofiler <- function(mutdata=x, esize=y, signatureset, kmin=2, kmax=4, ...){

	# load pathway dataset (already converted from .gmt file)
		edg <- signatureset

	# correct number of mutations for gene length (returns gene mutation scores)-without cagenes subset
		mutlength <- corgenelength_mod(x=mutdata, exonsize=esize)

	# transform mutlength
		mutlength <- t(mutlength)

	# correct for patient-specific mutation rate
		# calculate mutation rate
			patmutrate <- apply(mutlength, 2, sum)
			patmut0 <- which(patmutrate==0)
		# remove patients with mutationrate==0
			if(length(patmut0)>0){
				mutlength <- mutlength[,-patmut0,drop=F]
				patmutrate <- patmutrate[-patmut0]
			}

		# correct for mutation rate
			mutrate <- mutlength
		  for (p in 1:ncol(mutlength)){
		    mutrate[,p] <- mutlength[,p]/patmutrate[p]
		  }

  # correct gene scores for the number of pathways each gene belongs to
	  mutrate <- mutrate[which(row.names(mutrate) %in% colnames(edg)),]
	  genefreq <- apply(edg,2,sum)
	  genefreq <- genefreq[which(names(genefreq) %in% row.names(mutrate))]
	  mutratecor <- mutrate/genefreq

	# summarize gene mutation scores into pathway mutation scores
		signpat <- desparsify(edgx=edg, mutratecorx=mutratecor)

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
