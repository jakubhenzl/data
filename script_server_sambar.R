#oprava germinalnich dat podle checkgenesymbols
library(HGNChelper)
library(tidyverse)

germ_matrix<-read.csv2("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/sign_mat.csv",row.names = 1)
approved_genes<-checkGeneSymbols(colnames(germ_matrix))
colnames(germ_matrix)<-approved_genes$Suggested.Symbol
germ_matrix<-germ_matrix[,!is.na(colnames(germ_matrix))]


#vektor exon.size ze sambar
exon<-load("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/exon.size.RData")
approved_exon_size<-checkGeneSymbols(names(exon.size))


exon_size_correct<-exon.size
names(exon_size_correct)<-approved_exon_size$Suggested.Symbol
exon_size_correct<-exon_size_correct[!is.na(names(exon_size_correct))]

#nacteni drahove matice z gprofileru a kontrola genu pomoci checkgenesymbols
gprofiler_dataset<-read_csv2("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/gprofiler_dataset_ok.csv.gz")
gprofiler_dataset<-as.data.frame(gprofiler_dataset)
rownames(gprofiler_dataset)<-gprofiler_dataset$X1
gprofiler_dataset<-gprofiler_dataset[,-1]

approved_gprofiler<-checkGeneSymbols(colnames(gprofiler_dataset))

colnames(gprofiler_dataset)<-approved_gprofiler$Suggested.Symbol
gprofiler_dataset<-gprofiler_dataset[,!is.na(colnames(gprofiler_dataset))]


#sambar opravenej
source("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/corgenelength_mod.R") #corgenelength
source("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/sambar_gprofiler.R") #sambar_gprofiler
sambar_table<-sambar_gprofiler(mutdata=germ_matrix,esize=exon_size_correct,signatureset=gprofiler_dataset,kmin=2,kmax=7)
sambar_table<-as.data.frame(sambar_table)
write.csv(sambar_table,"/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/sambar_results")
