library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(Rgraphviz)
rm(list=ls())
library(grid)

M1=c("RPL11",	"MLLT11",	"DHCR24",	"ADIPOR1",	"ACTA1",	"TTC31",	"POLE4",	"PARD3B",	"GPC1",	"MPV17",	"ASB18",	"LRRC3B",	"SLC4A7",	"RHOA",	"KIAA1109",	"FBXL5",	"FLJ13197",	"PRR16",	"FLJ44896",	"FOXQ1",	"PPP1R11",	"PPIL1",	"YIPF3",	"HEY1",	"JAK2",	"TYRP1",	"CREB3",	"DCTN3",	"C10orf110",	"YME1L1",	"PSAP",	"CD151",	"HSD17B12",	"SNX15",	"ARHGAP42",	"NCAM1",	"ADAMTS15",	"PDDC1",	"TMEM9B",	"PATL1",	"RAB6A",	"SP1",	"PIP4K2C",	"TAS2R13",	"TAS2R19",	"PPP1R1A",	"LNX2",	"C14orf1",	"PLDN",	"SEPHS2",	"G6PC3",	"PFN1",	"CWC25",	"DCAKD",	"SOCS3",	"THOC4",	"STX10",	"NR2F6",	"RWDD2B",	"DGCR14",	"SF3A1",	"TMEM184B",	"PRAF2"
)
M2=c("C1orf144",	"IL23R",	"MAEL",	"RRAGC",	"ADIPOR1",	"SLC30A1",	"APEH",	"SLMAP",	"C3orf26",	"KBTBD12",	"TMEM115",	"FBXL5",	"ERGIC1",	"PPP2CA",	"PPP1R11",	"HIST1H2BK",	"CCDC25",	"AZIN1",	"TYRP1",	"NFX1",	"PDCL",	"ADO",	"YME1L1",	"DHX32",	"TMEM138",	"MED17",	"NCAM1",	"HPS5",	"UEVLD",	"PATL1",	"PIP4K2C",	"C12orf52",	"GTF2H3",	"SPRYD3",	"RB1",	"C14orf118",	"TTC8",	"TRMT61A",	"PLDN",	"SEPHS2",	"EXOC7",	"ZNF304",	"KCNE1",	"MALAT1"																				
)
M3=c("C1orf144",	"C1orf228",	"HECTD3",	"MEF2D",	"SMYD5",	"TMEM87B",	"PCNP",	"LMLN",	"NPRL2",	"TMEM115",	"PRR16",	"ERGIC1",	"PITX1",	"CAMK2A",	"TMEM151B",	"VPS52",	"HEY1",	"PSAP",	"C11orf80",	"GAB2",	"SPRYD3",	"RASSF9",	"CCDC92",	"FAM65A",	"VAC14",	"G6PC3",	"PFN1",	"TRAPPC1",	"VAT1",	"DCAKD",	"CD300LD",	"CHERP",	"TMEM161A",	"YIF1B",	"GSK3A",	"TTC28",	"SF3A1",	"PATZ1",	"PORCN",	"PRAF2",	"ZDHHC9"																							
)
M4=c("RNF207",	"CROCCP2",	"PARD3B",	"POMC",	"C3orf26",	"PRKCI",	"IP6K2",	"FLJ44896",	"SUMO4",	"GPLD1",	"SH2B2",	"DMRT3",	"DEAF1",	"PCSK7",	"CHD3",	"SOCS3",	"U2AF1"																																															
)
M5=c("RPL11",	"FOXD3",	"TMEM214",	"ASB18",	"RBM5",	"TTC14",	"B4GALT4",	"CDK14",	"CCDC25",	"TAS2R31",	"PEG3"																																																					
)
M6=c("HECTD3",	"C6orf208",	"SH2B2",	"C8orf56",	"PDCL",	"ACSF3",	"ZNF610"																																																									
)
M7=c("CCDC12",	"TBCEL",	"GAB2",	"PIP4K2C",	"APPL2",	"STX10",	"NR2F6",	"TBL1X"																																																								
)
M8=c("ASB18",	"TRPC1",	"CCDC12",	"CCAR1",	"HSD17B12",	"TBCEL",	"GAB2",	"PLDN",	"NR2F6"																																																							
)
M9=c("RNF207",	"RPL11",	"AKR7A2",	"TBCEL",	"SOCS3",	"C21orf82"																																																										
)
M10=c("CCNL1",	"LDHA",	"PPP1R1A",	"NDUFA13"																																																												
)
M11=c("SFRS18",	"BAALC",	"MLLT10",	"DDX17"																																																												
)
M12=c("C9orf130",	"MYST1"																																																														
)

M1=bitr(M1,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M1=bitr(M1,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M2=bitr(M2,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M3=bitr(M3,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M4=bitr(M4,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")

M5=bitr(M5,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")

M6=bitr(M6,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M7=bitr(M7,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M8=bitr(M8,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")

M9=bitr(M9,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")

M10=bitr(M10,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M11=bitr(M11,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")
M12=bitr(M12,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL","REFSEQ"),OrgDb = "org.Hs.eg.db")


list7<-list(M1=M1$ENTREZID,M2=M2$ENTREZID,M3=M3$ENTREZID,M4=M4$ENTREZID,M5=M5$ENTREZID,M6=M6$ENTREZID,M7=M7$ENTREZID,M8=M8$ENTREZID,M9=M9$ENTREZID,M10=M10$ENTREZID,M11=M11$ENTREZID)
# yy <- compareCluster(list7, fun="enrichKEGG", organism="hsa",pvalueCutoff=0.5)
xx <- compareCluster(list7,
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db", 
                     ont= "BP",
                     pvalueCutoff=0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.01)


yy <- compareCluster(list7, fun="enrichKEGG",organism="hsa", pvalueCutoff=0.05)

# as.data.frame(xx)
# plot(xx, type="dot", caption="KEGG Enrichment Comparison")
dotplot(xx,showCategory=10,includeAll=TRUE)
dotplot(yy,showCategory=10,includeAll=TRUE)
