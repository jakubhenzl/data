#oprava germinalnich dat podle checkgenesymbols
library(HGNChelper)
library(tidyverse)
library(SAMBAR)

germ_matrix<-read.csv2("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/sign_mat.csv",row.names = 1)
current_gene_map<-getCurrentHumanMap()
approved_genes<-checkGeneSymbols(colnames(germ_matrix),map=current_gene_map)
approved_genes[is.na(approved_genes$Suggested.Symbol), 3] <- approved_genes[is.na(approved_genes$Suggested.Symbol), 1]
colnames(germ_matrix)<-approved_genes$Suggested.Symbol

#vektor exon.size ze sambar

approved_exon_size<-checkGeneSymbols(names(SAMBAR::exon.size),map=current_gene_map)
approved_exon_size[is.na(approved_exon_size$Suggested.Symbol),3]<-approved_exon_size[is.na(approved_exon_size$Suggested.Symbol),1]
table(is.na(approved_exon_size$Suggested.Symbol))
exon_size_correct<-SAMBAR::exon.size
names(exon_size_correct)<-approved_exon_size$Suggested.Symbol


#nacteni drahove matice z gprofileru a kontrola genu pomoci checkgenesymbols
gprofiler_dataset<-read_csv2("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/gprofiler_dataset_ok.csv.gz")
gprofiler_dataset<-as.data.frame(gprofiler_dataset)
rownames(gprofiler_dataset)<-gprofiler_dataset$X1
gprofiler_dataset<-gprofiler_dataset[,-1]

approved_gprofiler<-checkGeneSymbols(colnames(gprofiler_dataset),map=current_gene_map)
approved_gprofiler[is.na(approved_gprofiler$Suggested.Symbol),3]<-approved_gprofiler[is.na(approved_gprofiler$Suggested.Symbol),1]

colnames(gprofiler_dataset)<-approved_gprofiler$Suggested.Symbol

#sambar opravenej
source("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/corgenelength_mod.R") #corgenelength
source("/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/sambar_gprofiler.R") #sambar_gprofiler
sambar_table<-sambar_gprofiler(mutdata=germ_matrix,esize=exon_size_correct,signatureset=gprofiler_dataset,kmin=2,kmax=7)
sambar_table<-as.data.frame(sambar_table)
write.csv(sambar_table,"/media/hdc/tmp/PCAWG_somatic_germline/Jakub/data/sambar_results")
