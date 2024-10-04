# gse ---> object with final results

r3 = c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2",
       "REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
       "REACTOME_PLATELET_AGGREGATION_PLUG_FORMATION" )

p1 = filter(gse$C5, ID == "GOBP_CYTOPLASMIC_TRANSLATION")


plot.combine.list = filter(gse$C2, ID %in% r3)
plot.combine.list@result = rbind(plot.combine.list@result, p1@result)
plot.combine.list@geneSets = c(plot.combine.list@geneSets, p1@geneSets)

test = gseaplot2(plot.combine.list, geneSetID = 1:4, rel_heights = c(1, 0.2, 0.2), base_size = 19, pvalue_table = T) 
test[[1]]$theme$legend.justification <- "left"
test[[1]]$theme$legend.position <- c(0.57,0.95)
test

sink(file = "ggplot.slots.log")                                            
str(test)
sink()


# ---------------- To match Jen's code -- and include only the significant results ------------------

# p1 = filter(gse$C5, ID == "GOBP_CYTOPLASMIC_TRANSLATION")

reactome.names = c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2")
gobp.names = c("GOBP_CYTOPLASMIC_TRANSLATION",
               # "GOBP_PLATELET_AGGREGATION",
               "GOBP_REGULATION_OF_PLATELET_ACTIVATION")
# C2 -- reactome
# C5 -- GOBP

reactome.Filtered.Names = filter(gse$C2, ID %in% reactome.names)
gobp.Filtered.Names = filter(gse$C5, ID %in% gobp.names)

reactome.Filtered.Names@result = rbind(reactome.Filtered.Names@result, gobp.Filtered.Names@result)
reactome.Filtered.Names@geneSets = c(reactome.Filtered.Names@geneSets, gobp.Filtered.Names@geneSets)
test = gseaplot2(reactome.Filtered.Names, geneSetID = 1:3, rel_heights = c(1, 0.2, 0.2), base_size = 19, pvalue_table = F) 
test[[1]]$theme$legend.justification <- "left"
test[[1]]$theme$legend.position <- c(0.57,0.95)

pdf(file = "sig-only-platelets-combined.pdf", width = 19, height = 19)
test
dev.off()

