
library(tidyverse)

germ_matrix<-read.csv("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/gprofiler_germ_despar_complete.csv",row.names = 1)
source("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/h_clustering_gprofiler.R")
results_2_15<-h_clustering_gprofiler(germ_matrix,kmin=2,kmax=15)
write.csv(results_2_15,"/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/results_2_15.csv")

germ_matrix<-t(germ_matrix)
x<-colSums(germ_matrix==0) 
y<-x[x>436]
not_want_pathways<-as.vector(names(y))
germ_filter<- germ_matrix %>% select(-not_want_pathways)
germ_filter<-t(germ_filter)

results_2_15_filtered<-h_clustering_gprofiler(germ_filter,kmin=2,kmax=15)
write.csv(results_2_15_filtered,"/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/results_2_15_filtered.csv")