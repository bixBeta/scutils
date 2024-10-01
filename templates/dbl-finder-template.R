# load in sobj ----

sobj.sample <- readRDS("${harmony_sobj}")

# create a v3 assay to convert to sce ----

sobj.sample[["RNA3"]] <- as(sobj.sample[["RNA"]], Class = "Assay")
DefaultAssay(sobj.sample) <- "RNA3"

# create sce from sobj.sample
sce = as.SingleCellExperiment(sobj.sample, assay = "RNA3")

# run scDblFinder ----
sce <- scDblFinder(sce, samples = "orig.ident", clusters = "harmony_clusters_res1")

table(sce$scDblFinder.class)

sce.meta = colData(sce) |> data.frame()

write.csv(x = sce.meta, file = "scDBL_doublet_identities.csv", quote = F)

cells2keep <- rownames(sce.meta |> filter(scDblFinder.class == "singlet"))

sobj.sample2 <- subset(sobj.sample, cells = cells2keep)

saveRDS(sobj.sample2, file = paste0("${id}", "_doubletFiltered_Integrated.RDS"))


