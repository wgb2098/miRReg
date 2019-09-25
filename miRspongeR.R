
library("miRspongeR")

ExpDatacsv  <- read.csv('miRNA&mRNA.csv', sep = ",", header = FALSE)#.csv the file is miRNA and mRNA exprofile

miR2Target  <- read.csv('miRNA_mRNA_matrix.csv', sep = ",", header = T)#.csv  this fiele is  miRNA and mRNA regulate matrix

pcceRInt <- spongeMethod(miR2Target, ExpDatacsv, method = "pc")
pcceRInt
#save the miRNA spongeT result.
write.csv( pcceRInt,"miRspongeRresult.csv",row.names=FALSE)
spongenetwork_Cluster <- netModule(pcceRInt[, 1:2])
sponge_Module_FEA <- moduleFEA(spongenetwork_Cluster)
Groundtruthcsv <- system.file("extdata", "Groundtruth.csv", package="miRspongeR")
Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",")
spongenetwork_validated <- spongeValidate(pcceRInt1[, 1:2], directed = FALSE, Groundtruth)

