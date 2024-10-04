library(dplyr)
library(purrr)
library(reshape2)
library(ggplot2)
nde = function(de.list, padj){
  
  
  z = list()
  
  for (i in 1:length(de.list)) {
    z[[i]] = length(which(pluck(de.list, i, "p_val_adj") < padj))
    names(z)[[i]] <- names(de.list)[[i]]
    
  }
  
  
  nDE.de.list = do.call("rbind",z)
  colnames(nDE.de.list) = deparse(substitute(de.list))
  
  
  return(nDE.de.list)
}


x.1 = nde(de.list = caseD1.controlD1.list, padj = 0.05)
x.2 = nde(de.list = caseD2.controlD2.list, padj = 0.05)
x.3 = nde(de.list = controlD1.controlD2.list, padj = 0.05)
x.4 = nde(de.list = caseD1.caseD2.list, padj = 0.05)

df = cbind(x.1,x.2,x.3,x.4)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(reshape2)
nde.gg = melt(df)
colnames(nde.gg) = c("cluster","condition", "nDE")

ggplot(nde.gg, aes(x=cluster, y=nDE , size = nDE , color = condition)) + scale_colour_viridis_d("condition") +
  facet_wrap(~condition) + ggtitle("Number of DE genes - SV4 - Tcells res 0.8") +
  geom_point(alpha=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

