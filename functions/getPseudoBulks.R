getByCluster = function(bulk_, cluster_){
  
  a = bulk_ |> select(matches(paste0("_", cluster_, "$"))) 
  return(a)
  
  
}

getPseudoBulk = function(sobj_, res_){
  
  bulk <- AggregateExpression(sobj_, group.by = c("orig.ident", paste0("harmony_clusters_res", res_)), return.seurat = F)
  bulk.matrix = as.data.frame(bulk$RNA) 
  
  nclusters = length(unlist(map(strsplit(colnames(bulk.matrix), split = "_"), ~pluck(.,2))) |> unique())

  pb.list = list()
  
  for (i in 1:nclusters) {
    message(paste("Getting Cluster:", i-1))
    pb.list[[i]] <-   getByCluster(bulk_ = bulk.matrix, cluster_ = i-1)
    names(pb.list)[[i]] <- paste0("cluster__", i-1, "__res", res_)
    
  }
  
  
  return(pb.list)
  
  
}


x = read.csv("sobjs.list.csv", header = F)

all.pbs = list()

for (i in 1:nrow(x)) {
  
  batch_ = (strsplit(x[i,1], split = "_")[[1]][2])
  message(paste0("Processing Batch: ", batch_ ))
  
  assign(x = paste0("sobj.batch.", batch_), value = readRDS(x[i,2]), envir = .GlobalEnv)
  
  all.pbs[[i]] <- getPseudoBulk(sobj_ = get(paste0("sobj.batch.", batch_)), res_ = 1.5)
  names(all.pbs)[[i]] <- paste0("batch_", batch_)
  
  rm(list = c(paste0("sobj.batch.", batch_)))

}

saveRDS(object = all.pbs, file = "Pseudobulk_RES.1.5.RDS")





