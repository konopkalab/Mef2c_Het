# Covariate plot
rm(list=ls())
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(biomaRt))

tab <- read.table("LM_EXON.txt")


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

MGI = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(tab) ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

comb <- merge(tab,MGI,by.x="row.names",by.y="MGI.symbol",all=T)
colnames(comb) <- c("MouseID","Estimate","Pval","FDR","Gene")

df <- comb %>% 
        mutate(Threshold = case_when(FDR < 0.05 ~ "DGE", FDR > 0.05  ~ "NOT_DGE")) %>%
        mutate(Direction = case_when(Estimate > 0 ~ "UpReg", Estimate < 0  ~ "DownReg")) 

# Load data for Data Base
psydge <- load(here("Human","PsychENCODE_DEGs.RData")) %>%
            get()
new_psydge_names <- names(psydge)
psydge <- map2(psydge, new_psydge_names, ~setnames(.x, 'Class', .y))

psymod <- load(here("Human","PsychEncode_Modules.RData")) %>%
            get()
new_psymod_names <- names(psymod)
psymod <- map2(psymod, new_psymod_names, ~setnames(.x, 'Mod', .y))

sfari <- load(here("Human","GeneSets_Disorders.RData")) %>%
            get()

# Mega list
psymod$CoolGenes=df

l <- c(psydge,psymod,sfari)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database=res[!(is.na(res$MouseID)),]


mef2cko <- load(here("Mouse","GeneSets_Mef2c_cKO.RData")) %>%
            get()
new_mef2cko<- names(mef2cko)
mef2cko <- map2(mef2cko, new_mef2cko, ~setnames(.x, 'Class', .y))

colnames(mef2cko[[1]])[1] <- "MouseID"
colnames(mef2cko[[2]])[1] <- "MouseID"
colnames(mef2cko[[3]])[1] <- "MouseID"


mef2cko$db <- database

l <- c(mef2cko)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="MouseID")
}, l)

database2=res[!(is.na(res$Estimate)),]

openxlsx::write.xlsx(database2, file = "Database.xlsx", colNames = TRUE, borders = "columns")






