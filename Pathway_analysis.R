.libPaths( c( .libPaths(), "/hpc/users/jayarp02/R/lib/4.0.3/") ) #  need to do this before loading any packages!!!

.libPaths()
library("ReactomePA")
load("RNASEQ_DEG.rData")
DEG_de = de[de$DEG=="DEG",]
subset = DEG_de
##############
#### Pathway analysis
##############

library("AnnotationDbi")
library("org.Hs.eg.db")
#columns(org.Hs.eg.db) # returns list of available keytypes
subset$entrez = mapIds(org.Hs.eg.db,
                       keys=subset$gene, #Column containing Ensembl gene ids
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals = "first")

write.csv(subset, file="AKI_COVID_DiffGenes_entrez.csv")

# ###POSLFC
# subset = DEG_de[DEG_de$LFC == "POSLFC",]
# 
# enriched_pathways<-enrichPathway(organism = "human",gene=unique(subset$entrez),pvalueCutoff=.1, readable=T,)
# x<-barplot(enriched_pathways, showCategory = 50, font.size = 20 )
# 
# descriptions <- as.character(x$data$Description)
# wrapped_text <-list()
# for (i in 1:length(descriptions)){
#   str1<- descriptions[i]
# 
#   wrapped_text[[i]]<-paste(strwrap(str1,width=60), collapse = "\n")
#   
#   
# }
# x$data$Description<- unlist(wrapped_text)
# x$data$Description<- factor(x$data$Description, levels = x$data$Description[order(-x$data$p.adjust)])
# x;
# 
# 
# ####NEGLFC
# subset = DEG_de[DEG_de$LFC == "NEGLFC",]
# 
# enriched_pathways<-enrichPathway(organism = "human",gene=unique(subset$entrez),pvalueCutoff=.1, readable=T,)
# x<-barplot(enriched_pathways, showCategory = 50, font.size = 20 )
# 
# descriptions <- as.character(x$data$Description)
# wrapped_text <-list()
# for (i in 1:length(descriptions)){
#   str1<- descriptions[i]
#   
#   wrapped_text[[i]]<-paste(strwrap(str1,width=60), collapse = "\n")
#   
#   
# }
# x$data$Description<- unlist(wrapped_text)
# x$data$Description<- factor(x$data$Description, levels = x$data$Description[order(-x$data$p.adjust)])
# x

### all LFC
enriched_pathways<-enrichPathway(organism = "human",gene=unique(subset$entrez),pvalueCutoff=.1, readable=T,)
x<-barplot(enriched_pathways, showCategory = 50, font.size = 20 )

descriptions <- as.character(x$data$Description)
wrapped_text <-list()
for (i in 1:length(descriptions)){
  str1<- descriptions[i]
  
  wrapped_text[[i]]<-paste(strwrap(str1,width=60), collapse = "\n")
  
  
}
x$data$Description<- unlist(wrapped_text)
x$data$Description<- factor(x$data$Description, levels = x$data$Description[order(-x$data$p.adjust)])
#x;




pdf(file = "pathways.pdf", width=60, onefile=T, height = 40)
##RESET
x;

### For GSEA need to get all proteins
geneList <- c(subset$logFC*-log10(subset$P.Value))
#geneList <- c((subset$FC))
names(geneList) <- subset$entrez
names(geneList)
geneList <- geneList[!names(geneList)==""]
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList, decreasing = T)
kk <- gsePathway(geneList ,pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 100 )
library(ggplot2)
ridge<- ridgeplot(kk, core_enrichment =T,label_format =  40 )
descriptions <- as.character(ridge$data$category)
wrapped_text <-list()
for (i in 1:length(descriptions)){
  str1<- descriptions[i]
  
  wrapped_text[[i]]<-paste(strwrap(str1,width=60), collapse = "\n")
  
  
}
ridge$data$category<- unlist(wrapped_text)
ridge$data$category<- factor(ridge$data$category, levels = unique(ridge$data$category[order(-ridge$data$value)]))
head(ridge$data)
ridge+ 
  xlab("Normalized Gene Expression")+
  ylab("Pathway")+
  labs(fill="FDR")+
  # scale_fill_brewer(palette="Set1")+
  theme(axis.title = element_text(size = 25),
        plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"),
        # plot.margin = margin(10, 10, 10, 70),
        axis.text.y = element_text(size = 25),
        axis.text.x= element_text(size = 15),
        
        # axis.text.x = element_text(angle = 45, hjust = 1, size=15),
        legend.text = element_text(size = 15),
        # legend.position = "none",
        title =element_text(size=20, face='bold'),
        # plot.title = element_text(size = 30,hjust = 0.5)
  )
dev.off()



