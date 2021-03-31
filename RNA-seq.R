rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
dat <- read.csv("gene_count_matrix.csv")
id <- dat[,1]
id <- as.data.frame(id)
eg <- bitr(id, fromType = "ENSEMBL",toType ="ENTREZID",OrgDb = "org.Mm.eg.db")
tab1 <- read.table("AB2.2-H4-1.tab",sep = "\t", header = T)
tab2 <- read.table("AB2.2-H4-2.tab",sep = "\t", header = T)
tab3 <- read.table("AB2.2-H4-3.tab",sep = "\t", header = T)
tab4 <- read.table("AB2.2-N4-1.tab",sep = "\t", header = T)
tab5 <- read.table("AB2.2-N4-2.tab",sep = "\t", header = T)
tab6 <- read.table("AB2.2-N4-3.tab",sep = "\t", header = T)
library(tidyverse)
id <- left_join(id,tab6, by=c("id"="Gene.ID"))
id <- id[,1:2]
id <- id[!duplicated(id),]
dat1 <- left_join(dat,id,by = c("gene_id"="id"))
dat1 <- dat1[,c(1,8,2:7)]
write.csv(dat1,file = "readcount_genename.csv")
count_tab <- read.csv("readcount_genename.csv", header = T)
count_tab = count_tab[,-c(1,3)]
row.names(count_tab) = count_tab[,1]
count_tab = count_tab[,-c(1)]
count_tab = count_tab[,c(4,5,6,1,2,3)]
colData <- data.frame(sample_id=colnames(count_tab), 
                      condition=factor(rep(c("Normal","Hypoxia"),each=3), levels=c("Normal","Hypoxia")))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_tab,
                              colData = colData,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="condition_Hypoxia_vs_Normal")
res <- res[order(res$padj), ]
resDF = as.data.frame(res)
resDF$gene_id = row.names(resDF)
resDF <- resDF[,c(7,1,2,3,4,5,6)]
DEG <- na.omit(resDF)
write.table(DEG, file = "DESeq2_DEG.tab", sep = "\t", quote = FALSE, row.names = FALSE)
genename <- read.csv("readcount_genename.csv", header = T)
genename <- genename[,2:3]
library(tidyverse)
DEG1 <- left_join(DEG,genename,by="gene_id")
save(colData, count_tab, dds, DEG, DEG1, genename, res, resDF, file = "HvsN.Rdata")
library(ggplot2)
library(ggrepel)
library(dplyr)
vol_plot <- function(DEG){
  logFC_cutoff <- 1
  DEG$change = as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                                ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                      '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                      '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
  )
  print(this_tile)
  print(head(DEG))
  library("ggplot2")
  g= ggplot(data=DEG, aes(x=log2FoldChange, y=-log10(padj), color=change)) +
    geom_point(alpha=0.4, size=3.5) +xlim(c(-6,6))+
    theme_set(theme_set(theme_bw(base_size=20)))+
    xlab("log2(fold change)") + ylab("-log10(p-value)") +
    geom_hline(yintercept = -log10(0.05), lty=2,color="black", lwd=0.8)+
    geom_vline(xintercept = c(-1,1),lty=2,color="black", lwd=0.8)+
    ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
    scale_colour_manual(values = c('blue','grey','red'))+theme_bw()+theme(plot.title = element_text(hjust=0.5),
                                                                          legend.position = "right",
                                                                          legend.title = element_blank())
   return(g)
}
pdf("volcano.pdf")
vol_plot(DEG1)
dev.off()
library(tidyverse)
library("pheatmap")
choose_gene <- head(rownames(DEG),50)
choose_matrix <- count_tab[choose_gene,]
choose_matrix$gene_id <- rownames(choose_matrix)
choose_matrix <- left_join(choose_matrix,genename,by="gene_id")
rownames(choose_matrix) <- choose_matrix$Gene.Name
choose_matrix <- choose_matrix[,-c(7,8)]
choose_matrix <- t(scale(t(choose_matrix)))
annotation_col <- data.frame(Condition=factor(rep(c("Normal","Hypoxia"),each=3),levels = c("Normal","Hypoxia")))
rownames(annotation_col) <- colnames(choose_matrix)
ann_colors = list(Condition = c(Normal = "light green", Hypoxia = "pink"))
pheatmap(choose_matrix, cluster_cols = F, annotation_col= annotation_col,
         cellwidth = 22, cellheight = 10,
         annotation_colors = ann_colors,
         filename = "heatmap.pdf")
library(ggplot2)
  library(DESeq2)
  vsd <- vst(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  PCAfig <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=4) +
    geom_text(aes(label=name), nudge_y = -1.0,nudge_x = -1.5)+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(PCAfig, filename = "HvsN_PCA.pdf")
library(tidyverse)
logFC_cutoff <- 1
DEG1$change = as.factor(ifelse(DEG1$padj < 0.05 & abs(DEG1$log2FoldChange) > logFC_cutoff,
                               ifelse(DEG1$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
Up_gene_number <- as.numeric(DEG1 %>% filter(change== "UP") %>% count())
Down_gene_number <- as.numeric(DEG1 %>% filter(change== "DOWN") %>% count())
Up_gene <- DEG1 %>% filter(change== "UP")
Down_gene <- DEG1 %>% filter(change== "DOWN")
write.csv(Up_gene,"Up_gene_list.csv")
save(DEG1,Up_gene,Down_gene, file = "HvsN_2.Rdata")
library(clusterProfiler)
library(org.Mm.eg.db)
gene_u <- bitr(Up_gene$gene_id,fromType = "ENSEMBL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)
ego <- enrichGO(gene= gene_u$ENTREZID,
                keyType = "ENTREZID",OrgDb = org.Mm.eg.db,
                ont= "ALL",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,
                readable      = TRUE)
write.csv(ego,file="GO_enrichment.csv")
save(ego,file = "ego_up.Rdata")
library(clusterProfiler)
library(org.Mm.eg.db) 
gene_d <- bitr(Down_gene$gene_id,fromType = "ENSEMBL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene= gene_d$ENTREZID,
                keyType = "ENTREZID",OrgDb = org.Mm.eg.db,
                ont= "ALL",pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,qvalueCutoff  = 0.05,
                readable      = TRUE)
write.csv(ego,file="GO_downregulated.csv")
save(ego,file = "ego_down.Rdata")
library(ggpubr)
library(ggthemes)
library(openxlsx)
library(ggsci)
data <- read.csv("GO_downregulated.csv")
if(T){
  library(tidyverse)
  data_BP <- data%>% filter(ONTOLOGY=="BP")%>% arrange(desc(Count))
  data_CC <- data%>% filter(ONTOLOGY=="CC")%>% arrange(desc(Count))
  data_MF <- data%>% filter(ONTOLOGY=="MF")%>% arrange(desc(Count))
  data_BP <- data_BP[1:10,]
  data_CC <- data_CC[1:10,]
  data_MF <- data_MF[1:10,]
  data <- na.omit(rbind(data_BP,data_CC, data_MF))
}
data$ONTOLOGY <- factor(data$ONTOLOGY, levels = c("MF","CC","BP"))
vline1 <- as.numeric(table(data$ONTOLOGY)[1])+0.5
vline2 <- as.numeric(table(data$ONTOLOGY)[1])+as.numeric(table(data$ONTOLOGY)[2])+0.5
data$Description<- sapply(as.character(data$Description),function(string) {ifelse (nchar(string)>40, paste(substr(string,1,40),"...",sep=""),string)})
a <- read.csv("GO_enrichment.csv")
data$Count <- -data$Count
data <- rbind(a,data)
data$ONTOLOGY <- factor(data$ONTOLOGY, levels = c("MF","CC","BP"))
vline1 <- as.numeric(table(data$ONTOLOGY)[1])+0.5
vline2 <- as.numeric(table(data$ONTOLOGY)[1])+as.numeric(table(data$ONTOLOGY)[2])+0.5
data$Description<- sapply(as.character(data$Description),function(string) {ifelse (nchar(string)>40, paste(substr(string,1,40),"...",sep=""),string)})
p1 <- ggbarplot(data,x="Description", y="Count", fill = "ONTOLOGY", rotate=T,
                sort.val = "asc", group="ONTOLOGY", sort.by.groups = T,
                xlab = "")+
  scale_fill_npg(alpha=0.8)+
  scale_y_continuous(expand = c(0,0),limits = c(-15,12))+geom_vline(xintercept = c(vline1,vline2),linetype="dashed", size=1)+
  geom_hline(yintercept = 0,size=1)+theme_base()
ggsave("GO_MERGE_Bar.pdf",p1, width = 12,height = 12, units = "in")
load("HvsN.Rdata")
genename <- read.csv("readcount_genename.csv", header = T)
meso <- read.table("meso_gene_LM.txt")
meso <- unique(meso)
colnames(meso) <- c("Gene.Name")
library(tidyverse)
meso_count <- left_join(meso,genename,by="Gene.Name")
rownames(meso_count) <- meso_count[,1]
meso_count <- meso_count[,-c(1:3)]
meso_count <- t(scale(t(meso_count)))
meso_count <- meso_count[,c(4:6,1:3)]
meso_count <- ifelse(meso_count>=-1.2&meso_count<=1.75, meso_count, ifelse(meso_count>1.75, 1.75, -1.2))
annotation_col <- data.frame(Condition=factor(rep(c("Normal","Hypoxia"),each=3),levels = c("Normal","Hypoxia")))
rownames(annotation_col) <- colnames(meso_count)
ann_colors = list(Condition = c(Normal = "light green", Hypoxia = "pink"))
library(pheatmap)
pheatmap(meso_count, cluster_cols = F, annotation_col= annotation_col,
         cellwidth = 22, cellheight = 10,
         annotation_colors = ann_colors,
         show_colnames = F,
         main = "Mesoderm",
         filename = "meso_gene.pdf")
endo <- read.table("endo_gene_LM.txt")
endo <- unique(endo)
colnames(endo) <- c("Gene.Name")
library(tidyverse)
endo_count <- left_join(endo,genename,by="Gene.Name")
rownames(endo_count) <- endo_count[,1]
endo_count <- endo_count[,-c(1:3)]
endo_count <- t(scale(t(endo_count)))
endo_count <- endo_count[,c(4:6,1:3)]
endo_count <- ifelse(endo_count>=-1.2&endo_count<=1.75, endo_count, ifelse(endo_count>1.75, 1.75, -1.2))
annotation_col <- data.frame(Condition=factor(rep(c("Normal","Hypoxia"),each=3),levels = c("Normal","Hypoxia")))
rownames(annotation_col) <- colnames(endo_count)
ann_colors = list(Condition = c(Normal = "light green", Hypoxia = "pink"))
library(pheatmap)
pheatmap(endo_count, cluster_cols = F, annotation_col= annotation_col,
         cellwidth = 22, cellheight = 10,
         annotation_colors = ann_colors,
         show_colnames = F,
         main = "Endoderm",
         filename = "endo_gene.pdf")
ecto <- read.table("ecto_gene_LM.txt")
ecto <- unique(ecto)
colnames(ecto) <- c("Gene.Name")
library(tidyverse)
ecto_count <- left_join(ecto,genename,by="Gene.Name")
rownames(ecto_count) <- ecto_count[,1]
ecto_count <- ecto_count[,-c(1:3)]
ecto_count <- t(scale(t(ecto_count)))
ecto_count <- ecto_count[,c(4:6,1:3)]
ecto_count <- ifelse(ecto_count>=-1.2&ecto_count<=1.75, ecto_count, ifelse(ecto_count>1.75, 1.75, -1.2))
annotation_col <- data.frame(Condition=factor(rep(c("Normal","Hypoxia"),each=3),levels = c("Normal","Hypoxia")))
rownames(annotation_col) <- colnames(ecto_count)
ann_colors = list(Condition = c(Normal = "light green", Hypoxia = "pink"))
library(pheatmap)
pheatmap(ecto_count, cluster_cols = F, annotation_col= annotation_col,
         cellwidth = 22, cellheight = 10,
         annotation_colors = ann_colors,
         show_colnames = F,
         main = "Ectoderm",
         filename = "ecto_gene.pdf")
load("HvsN.Rdata")
DEG2 <- DEG1[order(-DEG1$log2FoldChange),]
head(DEG2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(enrichplot)
library(ReactomePA)
eg <- bitr(DEG2$gene_id,
           fromType = "ENSEMBL",
           toType = "ENTREZID",
           OrgDb = "org.Mm.eg.db")
head(eg)
DEG2 <- inner_join(eg,DEG2,by=c("ENSEMBL"="gene_id"))
head(DEG2)
genelist <- DEG2$log2FoldChange
names(genelist) <- DEG2$ENTREZID
kk <- gseKEGG(geneList = genelist,
              organism = "mmu",
              nPerm = 1000,
              minGSSize = 10,
              maxGSSize = 1000,
              pvalueCutoff = 1,
              verbose = FALSE)
go <- gseGO(genelist, 'org.Mm.eg.db', keyType = "ENTREZID", ont="BP", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
gsea_plot <- function(go_item){
  gseaplot(go,geneSetID=go$ID[str_detect(go$Description,paste("^",go_item,"$",sep = ""))],
           title=go$Description[str_detect(go$Description,paste("^",go_item,"$",sep = ""))])
  ggsave(paste(go_item,".pdf",sep = ""),plot = last_plot(),width = 8, height = 8)
}
gsea_plot("regulation of non-canonical Wnt signaling pathway")
gseaplot(kk,geneSetID = "mmu04066",title = kk$Description[str_detect(kk$Description,"^HIF-1 signaling pathway$")])
ggsave("HIF-1 signaling pathway.pdf",plot = last_plot(),width = 12, height = 12)
gseaplot(kk,geneSetID = "mmu04310",title = kk$Description[str_detect(kk$Description,"^Wnt signaling pathway$")])
ggsave("Wnt signaling pathway.pdf",plot = last_plot(),width = 12, height = 12)



