corgenelength_mod <- function(x, exonsize, ...){

	# match with mutation data
		x <- x[,order(colnames(x))]
		x <- x[,which(colnames(x) %in% names(exonsize))]
		exonsize <- exonsize[which(names(exonsize) %in% colnames(x))]

	# normalize the number of mutations by gene length
		y <- sweep(x,2,exonsize,'/')

	return(y) # y is mutrate
}
