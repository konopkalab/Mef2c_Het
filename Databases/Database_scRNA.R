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

df <- tab %>% 
		rownames_to_column("Gene") %>%
        mutate(Threshold = case_when(FDR < 0.05 ~ "DGE", FDR > 0.05  ~ "NOT_DGE")) %>%
        mutate(Direction = case_when(Estimate > 0 ~ "UpReg", Estimate < 0  ~ "DownReg")) 


scRNA <- load(here("Mouse","GeneSets_scMouse.RData")) %>%
            get()
new_mef2cko<- names(scRNA)
scRNA <- map2(scRNA, new_mef2cko, ~setnames(.x, 'Class', .y))

microglia <- load(here("Mouse","GeneSets_Microglia.RData")) %>%
            get()
new_mef2cko<- names(microglia)
microglia <- map2(microglia, new_mef2cko, ~setnames(.x, 'Class', .y))

microglia2 <- load(here("Mouse","Neuron_Micro.RData")) %>%
            get()
new_mef2cko<- names(microglia2)
microglia2 <- map2(microglia2, new_mef2cko, ~setnames(.x, 'Class', .y))

scRNA$db <- df

l <- c(scRNA,microglia,microglia2)
res <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, l)

database2=res[!(is.na(res$Estimate)),]

openxlsx::write.xlsx(database2, file = "Database_scRNA.xlsx", colNames = TRUE, borders = "columns")






