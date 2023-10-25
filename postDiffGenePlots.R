#!/usr/bin/env Rscript

.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") ) #  need to do this before loading any packages!!!
.libPaths()
library("BiocManager")
#BiocManager::install("ReactomePA", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#BiocManager::install("ggridges", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
#BiocManager::install("ggupset", lib="/hpc/users/jayarp02/R/lib/4.0.3/")
library("ReactomePA")
library("ggridges")
library("gdata")
library("ggupset")
library("tidyverse")
library("AnnotationDbi")
library("org.Hs.eg.db")
#library("ensembldb")
#library("EnsDb.Hsapiens.v86")
library("biomaRt")
library("data.table")
library("dplyr")
library("ggplot2")
library("httr")

new_config <- httr::config(ssl_verifypeer = FALSE)
verify_peer <- TRUE
httr::set_config(new_config, override = FALSE)

options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)

file = paste0(args[1], "/RNASEQ_DEG.rData")
print(file)

load(file)

de = de[adj.P.Val<0.05]
colnames(de)
ttable = de

write.csv(file=paste0(args[1],"/RNASEQ_DEG.csv"), de)

### Map gene symbols
#mapping = read.csv(file="ensembl_entrez_mapping.txt", stringsAsFactors = F)
#head(mapping)

### map ensembl ids to gene symbols
#ttable$gene_name <- mapping$NCBI.gene..formerly.Entrezgene..ID[match(ttable$gene, mapping$Gene.stable.ID)]
#ttable$gene_symbol <- mapping$Gene.stable.ID[match(ttable$gene, mapping$Gene.stable.ID)]
#ttable$gene_description <- mapping$NCBI.gene..formerly.Entrezgene..description[match(ttable$gene, mapping$Gene.stable.ID)]

### ensembl gene IDs to HGNC gene symbols
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org" )
mapped_genes <- data.frame(getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                                 filters = "ensembl_gene_id", values = as.vector(ttable$gene),
                                 mart = mart))



#listAttributes(mart)
mapped_genes

ttable$gene_name <- mapped_genes$entrezgene_id[match(ttable$gene, mapped_genes$ensembl_gene_id)]
saveRDS(ttable, file="enriched_genes.rds")
geneList <- c(ttable$logFC)
names(geneList) <- ttable$gene_name
geneList <- geneList[!names(geneList)==""]
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList, decreasing = T)

pdf(paste0(args[1],"/diff_geneexpr_enrichment_analysis.pdf"), height=30, width=25)
print(geneList)

kk <- gsePathway(geneList ,pvalueCutoff = 0.08, minGSSize = 10  )

ridge<- ridgeplot(kk, core_enrichment = T)

ridge+ 
  # scale_fill_brewer(palette="Set1")+
  theme(axis.title = element_text(size = 25),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        # plot.margin = margin(10, 10, 10, 70),
        # axis.text.y = element_text(size = 20),
        axis.text.x= element_text(size = 15),
        
        # axis.text.x = element_text(angle = 45, hjust = 1, size=15),
        legend.text = element_text(size = 20),
        # legend.position = "none",
        title =element_text(size=20, face='bold'),
        # plot.title = element_text(size = 30,hjust = 0.5), 
        plot.background = element_rect(fill="transparent"))

print(kk$Description)

library(enrichplot)
gseaplot(kk, geneSetID = 6, by = "runningScore", title = kk$Description[3],
)
upsetplot(kk)
gseaplot2(kk, geneSetID = 6, title = kk$Description[3],
)



##############
#### Pathway analysis
##############
#ttable$entrez <- protein_metadata$EntrezGeneID[match(rownames(ttable), protein_metadata$TargetFullName)]
subset <- ttable[ttable$logFC>0,]
enriched_pathways<-enrichPathway(gene=subset$gene_name,pvalueCutoff=0.1, readable=T)
x<-barplot(enriched_pathways, showCategory = 15, font.size = 15, )
# write.csv(x$data, file ="./pathways/PTL_CD8_N_Cnt_CD8_N_upregulated_enriched_reactome_pathways.csv" )
descriptions <- as.character(x$data$Description)
wrapped_text <-list()
for (i in 1:length(descriptions)){
  str1<- descriptions[i]
  # wrapped_text[[i]]<-paste(sapply(seq(1, nchar(str1), 16), function(i) paste0(substring(str1, i, min(i + 15, nchar(str1))), '\n')), collapse='')
  wrapped_text[[i]]<-paste(strwrap(str1,width=60), collapse = "\n")
  
  
}
x$data$Description<- unlist(wrapped_text)
x$data$Description<- factor(x$data$Description, levels = x$data$Description[order(-x$data$p.adjust)])
x
dev.off()
