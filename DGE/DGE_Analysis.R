# Load libraries
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(preprocessCore))

####################################
# Load Inputs and filter for Neun+ #
####################################
RPKM=read.table("RPKM_EXON.txt",sep="\t",header=T)
COUNT=read.table("COUNT_EXON.txt",sep="\t",header=T)

first=apply(RPKM, 1, function(x) (all(x[1:4] >= 0.5) | all(x[5:8] >= 0.5)))
df <- COUNT[first,]

##########################
#    normalization       #
##########################
dat <- log2(cpm(df)+1)

##########################
# Quantile normalization #
##########################
p <- normalize.quantiles(as.matrix(dat))
rownames(p) <- rownames(dat)
colnames(p) <- colnames(dat)
#write.table(p, "OUTPUTS_GENEBODY_NEUN/NeuN_Primates_GeneBody_CPM.txt",sep="\t",quote=F)

pd=data.frame(row.names = colnames(p), Treatment=c(rep("HET",4),rep("WT",4)))


#########################
#         PCA           #
#########################
pdf("PCA.pdf",width=5,height=5,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x,Treatment=pd$Treatment)
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",color = "Treatment",palette=c("steelblue","darkgrey","green"), shape = 21, size = 3)+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic() 
dev.off()

########################
# Human vs Chimpanzee  #
########################
TRAITSfilt <- droplevels(pd)
Data=t(p)
output <- data.frame(matrix(nrow=ncol(Data), ncol=3, dimnames = list(colnames(Data), c("Estimate", "Pval", "Warning"))))
output[,] <- NA
for (i in 1:ncol(Data)) {   
            Model=tryCatch(lm(as.formula(paste("Data[,i] ~ ", paste(colnames(TRAITSfilt),collapse = " + "))), data = TRAITSfilt),warning =  function(w) w)
                if (i %% 1000 == 0) {cat(paste("Done on gene ",i,"\n",sep=""))}
                if(typeof(Model) == "list"){
                    coefs = data.frame(coef(summary(Model)))
                    t_value = coefs["Treatment", "t.value"]
                    output[i,"Pval"] = 2 * (1 - pnorm(abs(t_value)))
                    output[i,"Estimate"]= -1 * coefs["Treatment", "Estimate"]
                    } else {
                    output[i,"Warning"] = as.character(Model)
                    output[i, "Estimate"] = 0
                    output[i,"Pval"] = 1
            }
        }

DGE <- output
DGE$Warning <- NULL
DGE$FDR <- p.adjust(DGE$Pval,"BH")
write.table(DGE, "LM_EXON.txt",sep="\t",quote=F)
xlsx::write.xlsx(DGE, file="LM_EXON.xlsx",sheetName = "EXON DEG",row.names=TRUE, showNA=FALSE)


sign <- DGE[DGE$FDR < 0.05,]

mat <- p[rownames(p)%in% rownames(sign),]
colnames(mat) <- c("Het_1","Het_2","Het_3","Het_4","Wt_1","Wt_2","Wt_3","Wt_4")
anno <- data.frame(row.names = colnames(mat), Treatment=c(rep("HET",4),rep("WT",4)))
Treatment        <- c("red", "black")
names(Treatment) <- c("HET", "WT")
anno_colors <- list(Treatment = Treatment)
pdf("Heatmap.pdf",width=4,height=6)
pheatmap(mat,scale="row",show_rownames = F,annotation=anno,annotation_colors = anno_colors)
dev.off()

# Convert to human 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(sign) ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

signComb <- merge(sign,MGI,by.x="row.names",by.y="MGI.symbol",all=F)

df <- signComb %>%
                mutate(Direction = case_when(Estimate > 0 ~ "Upreg", Estimate < 0  ~ "Downreg")) %>%
                dplyr::select(HGNC.symbol,Direction) %>%
                dplyr::rename(Gene = HGNC.symbol)

tmp <- data.frame(Gene = df$Gene, Direction = rep("All",nrow(df)))

df <- rbind(df,tmp)

write.table(df,"DEGs_For_Enrichment.txt",sep="\t",quote=F,row.names=F)

# Vulcano Plot
tab <- read.table("LM_EXON.txt")
df <- tab %>% 
        rownames_to_column("Names") %>%
        mutate(Threshold = case_when(FDR < 0.05 ~ "TRUE", FDR > 0.05  ~ "FALSE")) %>%
        mutate(Direction = case_when(Estimate > 0 ~ "UpReg", Estimate < 0  ~ "DownReg")) %>%
        mutate(LOG = -log10(FDR), ABS = abs(Estimate)) 

df$LOG[!is.finite(df$LOG)] <- 12


top_labelled <- tbl_df(df) %>% 
                  group_by(Direction) %>% 
                  top_n(n = 10, wt = LOG)

pdf("Vulcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df, 
            x = "Estimate", 
            y = "LOG",
            color = "Threshold",
            palette=c("grey","red"),
            size = 1,
            alpha=0.3,
            shape=19)+
      xlab("log2(Fold Change)")+ 
      ylab("FDR")+
      geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) + 
      geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
      geom_text_repel(data = top_labelled, 
                      mapping = aes(label = Names), 
                      size = 5,
                      box.padding = unit(0.4, "lines"),
                      point.padding = unit(0.4, "lines"))+
      theme(legend.position="none")+
      ylim(0,20)+
      xlim(-2,2)
dev.off()













