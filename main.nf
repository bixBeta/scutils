nextflow.enable.dsl=2

params.sheet            = "sample-sheet.csv"
ch_sheet = channel.fromPath(params.sheet)


meta_ch = ch_sheet
                |  splitCsv( header:true )
                |  map { row -> [row.label, [file("RDS/" + row.sobj)]] }
                |  view



process HARMONY {

    label 'process_high'
    tag "$id"
    maxForks 3

    publishDir "FIGURES/${id}"      , mode: 'symlink', pattern: "*.png", overwrite: true
    publishDir "SOBJ_METAS/${id}"   , mode: 'symlink', pattern: "*.csv", overwrite: true
    publishDir "SOBJ_RDS/${id}"     , mode: 'symlink', pattern: "*.RDS", overwrite: true


    input:
        tuple val(id), path(sobj)

    output:
        path "*.png"
        path "*.csv"
        tuple val(id), path("*.RDS")                    , emit: harmony_sobj


    script:

    """
        #!/usr/bin/env Rscript
        .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2")

        library(Seurat)
        library(SeuratWrappers)
        library(ggplot2)
        library(patchwork)
        library(reticulate)
        library(glmGamPoi)
        library(dplyr)
        library(plotly)

        sobj <- readRDS("${sobj}")
        sobj.meta <- sobj@meta.data

        DefaultAssay(sobj) <- "SCT"

        var.feats = VariableFeatures(sobj)
        tcr.feats = var.feats[grep("^TR[A|D|B]", x = var.feats)]

        VariableFeatures(sobj) <- setdiff(var.feats, tcr.feats)
        sobj <- RunPCA(sobj, npcs = 50, verbose = T)

        png(paste0("${id}", "_ElbowPlot.png"), width = 1080, height = 1080)
        ElbowPlot(sobj, ndims = 50, reduction = "pca")
        
        sobj <- IntegrateLayers(object = sobj, method = HarmonyIntegration, orig.reduction = "pca", 
                                new.reduction = "harmony", normalization.method = "SCT",
                                verbose = T)

        sobj <- FindNeighbors(sobj, dims = 1:50, reduction = "harmony")
        
        sobj <- FindClusters(sobj, resolution = c(0.8,1,1.5,2), cluster.name = c("harmony_clusters_res0.8" , "harmony_clusters_res1", "harmony_clusters_res1.5", "harmony_clusters_res2"))
        
        sobj.meta = sobj@meta.data
        write.csv(sobj.meta, file = paste0("$id", "_metadata.csv"), quote = F)

        sobj <- RunUMAP(sobj, dims = 1:50, reduction = "harmony", reduction.name = "umap.harmony")

        Idents(sobj) <- sobj\$harmony_clusters_res0.8
        a = DimPlot(sobj, label = T, label.size = 7, raster = F, reduction = "umap.harmony", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 0.8")
        a2 = DimPlot(sobj, raster = F, reduction = "umap.harmony", pt.size = 0.3, group.by = "orig.ident") + ggtitle("Harmony Cluster Resolution : 0.8")

        Idents(sobj) <- sobj\$harmony_clusters_res1
        b = DimPlot(sobj, label = T, label.size = 5, raster = F, reduction = "umap.harmony", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 1")
        b2 = DimPlot(sobj,  raster = F, reduction = "umap.harmony", group.by = "orig.ident", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 1")
        
        Idents(sobj) <- sobj\$harmony_clusters_res1.5
        c = DimPlot(sobj, label = T, label.size = 5, raster = F, reduction = "umap.harmony", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 1.5")
        c2 = DimPlot(sobj, raster = F, reduction = "umap.harmony", group.by = "orig.ident", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 1.5")

        Idents(sobj) <- sobj\$harmony_clusters_res2
        d = DimPlot(sobj, label = T, label.size = 5, raster = F, reduction = "umap.harmony", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 2")
        d2 = DimPlot(sobj,  raster = F, reduction = "umap.harmony", group.by = "orig.ident", pt.size = 0.3) + ggtitle("Harmony Cluster Resolution : 2")

        png(paste0("${id}", "_Harmony-by-Sample-UMAP1-res0.8.png"), width = 2400, height = 900, res = 150)
        print( a | a2 )
        dev.off()

        png(paste0("${id}", "_Harmony-by-Sample-UMAP1-res1.png"), width = 2400, height = 900, res = 150)
        print( b | b2 )
        dev.off()

        png(paste0("${id}", "_Harmony-by-Sample-UMAP1-res1.5.png"), width = 2400, height = 900, res = 150)
        print( c | c2)
        dev.off()

        png(paste0("${id}", "_Harmony-by-Sample-UMAP2-res2.png"), width = 2400, height = 900, res = 150)
        print( d | d2 )
        dev.off()

        saveRDS(sobj, file = paste0("${id}", "_Harmony-by-sample-sobj.RDS"))

    """





}




process LOUPER {

 label 'process_loupe'
    tag "$id"
    maxForks 3
    publishDir "SOBJ_LOUPEr/${id}"     , mode: 'symlink', pattern: "*.cloupe", overwrite: true


    input:
        tuple val(id), path(harmony_sobj)

    output:
        path "*.cloupe"

    script:
    
    """
    #!/usr/bin/env Rscript
    .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2")    

    library(Seurat)
    library(SeuratWrappers)
    library(ggplot2)
    library(patchwork)
    library(reticulate)
    library(glmGamPoi)
    library(dplyr)
    library(plotly)
   
    loupeR::setup("yes")

    # import the library
    library("loupeR")

    # Gene Expression RNA assay, if working with SCT, please change the "RNA" to "SCT"
    sobj <- readRDS("${harmony_sobj}")

    assay <- sobj[["SCT"]]

    # get counts matrix from either the old or newer formats of assay
    counts <- counts_matrix_from_assay(assay)


    # convert the count matrix, clusters, and projections into a Loupe file
    create_loupe(
    counts,
    clusters = select_clusters(seurat_obj),
    projections = select_projections(seurat_obj),
    output_name = "${id}_sct_sobj"
    )


    """
}



process ADT{

    label 'process_high'
    tag "$id"
    maxForks 3
    publishDir "SOBJ_ADT/${id}"     , mode: 'symlink', pattern: "*.RDS", overwrite: true


    input:
        tuple val(id), path(harmony_sobj)

    output:
        path "*.RDS"

    script:
    
    """
    #!/usr/bin/env Rscript
    .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2")    


    library(Seurat)
    library(SeuratWrappers)
    library(ggplot2)
    library(patchwork)
    library(reticulate)
    library(glmGamPoi)
    library(dplyr)
    library(plotly)

    # Gene Expression RNA assay, if working with SCT, please change the "RNA" to "SCT"
    sobj.sample <- readRDS("${harmony_sobj}")

    DefaultAssay(sobj.sample) <- "ADT"

    sobj.sample[["ADT"]] <- JoinLayers(sobj.sample[["ADT"]])

    sobj.sample <- NormalizeData(sobj.sample, normalization.method = "CLR", margin = 2, assay = "ADT")


    saveRDS(sobj.sample, file = paste0("${id}", "_Harmony-by-sample-ADT-normalized-sobj.RDS"))





    """




}



process DBL{

    label 'process_scdbl'
    tag "$id"
    maxForks 3
    publishDir "SOBJ_DOUBLET_FILTERED_SOBJS/${id}"     , mode: 'symlink', pattern: "*.RDS", overwrite: true
    publishDir "SOBJ_DOUBLET_FILTERED_METADATA/"       , mode: 'symlink', pattern: "*.csv", overwrite: true

    input:
        tuple val(id), path(harmony_sobj)

    output:
        path "*.RDS"
        path "*.csv"                     , emit: scdbl_meta

    script:

        """

            #!/usr/bin/env Rscript
            .libPaths("/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2")    


            library(Seurat)
            library(SeuratWrappers)
            library(ggplot2)
            library(patchwork)
            library(reticulate)
            library(glmGamPoi)
            library(dplyr)
            library(plotly)            
            library(scDblFinder)


            # load in sobj ----

             sobj.sample <- readRDS("${harmony_sobj}")

            # create a v3 assay to convert to sce ----
    
             sobj.sample[["RNA3"]] <- as(sobj.sample[["RNA"]], Class = "Assay")
             DefaultAssay(sobj.sample) <- "RNA3"

            # create sce from sobj.sample
            sce = as.SingleCellExperiment(sobj.sample, assay = "RNA3")

            # run scDblFinder ----
            sce <- scDblFinder(sce, samples = "orig.ident", clusters = "harmony_clusters_res1")

            table(sce\$scDblFinder.class)

            sce.meta = colData(sce) |> data.frame()
            write.csv(x = sce.meta, file = paste0("${id}", "_scDBL_doublet_identities.csv"), quote = F)

            cells2keep <- rownames(sce.meta |> filter(scDblFinder.class == "singlet"))
            sobj.sample2 <- subset(sobj.sample, cells = cells2keep)

            saveRDS(sobj.sample2, file = paste0("${id}", "_doubletFiltered_Integrated.RDS"))
            saveRDS(sce, file = paste0("${id}",  "_doubletCalled_SCE_CLASS_Integrated.RDS"))



        """




}

workflow {
  
    HARMONY(meta_ch)


    harm_ch = HARMONY.out.harmony_sobj
            | view

    ADT(harm_ch)

    DBL(harm_ch)
}