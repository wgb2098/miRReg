library(GEOquery)
GPL5175 <-getGEO('GPL5175',destdir =".")
GPL5175_anno <-  data.table::fread('./GPL5175.soft¡ä,skip ="ID")
colnames(as.data.frame.table(GPL5175_anno))
ID2symbol=data.frame(GPL5175_anno)[,c(1£¬10)]
ID2symbol
library(dplyr)
library(tidyr)
probe2symbol_df <- GPL5175_anno
dplyr::select(probe2symbol_df,ID,gene_assignment)
dplyr::filter(probe2symbol_df,gene_assignment!="---")
tidyr::separate(probe2symbol_df,gene_assignment,c("drop","symbol"),sep="//")
dplyr::select(-drop)
ID3symbol=tidyr::separate(probe2symbol_df,gene_assignment,c("drop","symbol"),sep="//")
ID4symbol=dplyr::filter(ID3symbol,ID3symbol$symbol=="NA")
ID4symbol=dplyr::select(ID3symbol,ID,symbol)
ID4symbol