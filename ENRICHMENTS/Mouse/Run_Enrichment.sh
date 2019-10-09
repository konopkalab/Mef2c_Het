# Enrichment.r
# -g = Two columns table of input genes with specific association from your study
# -l = list of two columns tables with gene - disease association. E.g. Gene1 SYN
# -p = make a bubble chart with OR and -log10(FDR)
# -b = background (protein coding = 19776)
# -o = output label for statistics and viz
# -W/-H = width/height of the plot. 
Rscript Enrichment.r -g DEGs_For_Enrichment_Mouse.txt -l GeneSets_scMouse.RData -p -b 19776 -o scMouse -W 5 -H 5