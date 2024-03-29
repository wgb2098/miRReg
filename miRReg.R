####################################################################
# STEP 1 : miRNA-mRNA regulatory network and module 
####################################################################
# DE mRNA and DE miRNA 
####################################################################
rm(list=ls())
library(limma)
getwd()
exp<-read.table("gse35776.txt", header=T,sep="\t")
names(exp)
d <- exp[, 2:25]
rownames(d) <- exp[, 1]
# Pheno data file
pheno <- read.table("group2.txt",header = F)
colnames(pheno) <- c("sample","group")

Group<-factor(pheno$group,levels=levels(pheno$group))

design<-model.matrix(~0+Group)
y <- voom(d,design,plot=TRUE)
colnames(design)
fit <-lmFit(y,design)
##Designing Contrast Matrix for group Differentiation
cont.wt<-makeContrasts(Groupcontrol-Grouptof,levels=design)
fit2 <-contrasts.fit(fit,cont.wt)
fit3<-eBayes(fit2)
DE<-topTable(fit3,number = 'all')
write.csv(DE, "DE_mRNA.csv")
selected = rownames(DE[DE$P.Value <0.05,])
esetSel = y[selected, ]
heatmap(esetSel@.Data[[1]])
###############################################################
# LASSO regession algorithm
###############################################################
datafile <- read.csv("mrn.csv")
datafile<-t(datafile)
write.table(datafile,"mrn-1.csv",sep=",", row.names=FALSE)
datafile1 <- read.csv("mirna.csv")
datafile1<-t(datafile1)
write.table(datafile1,"mirna-1.csv",sep=",", row.names=FALSE)

datafile <- read.csv("mrn-new.csv")
datafile1 <- read.csv("mirna-new.csv")

x <- as.matrix(datafile1[, 1:149])  # Convert to matrix

for ( i in 1:ncol(datafile))

{ y <- as.matrix(datafile[, i ])
  gla <- cv.glmnet(x, y, nfolds =5)
  gla.best <- gla$glmnet.fit  # Corresponding best model
  gla.coef <- coef(gla$glmnet.fit, s = gla$lambda.1se)  # Extract model coefficients
  savename=paste("A",i,".csv",sep="")
  write.csv(as.matrix(t(abs(gla.coef))),savename,row.names = FALSE)
}
############# synthesis files
files.name=list.files(pattern = "csv");
files.length=length(files.name);
newdata=read.csv('A1.csv',head=T,sep=",");
newdata<-as.matrix(newdata);
for(i in 2:files.length)
{
  outname=paste("A",i,".csv",sep="")
  tmp=read.csv(outname,head=T,sep=",");
  tmp=as.matrix(tmp);
  newdata=rbind(newdata, tmp[1,]);
 }
# save the results
write.csv(newdata,"result.csv",row.names = FALSE)
#########################################################################
#id  change gene symbol
########################################################################
library(GEOquery)
GPL5175 <-getGEO('GPL5175',destdir =".")
GPL5175_anno <-  data.table::fread('./GPL5175.soft��,skip ="ID")
colnames(as.data.frame.table(GPL5175_anno))
ID2symbol=data.frame(GPL5175_anno)[,c(1��10)]
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
#########################################################################
# biclique miRNA-mRNA
##########################################################################
rm(list=ls()) 
library(igraph)
library(parallel)
library(Matrix)
library(glmnet)
library(plyr)
dircode   = "E:\\R\\okok\\sourcecode"
dirdata   = "E:\\R\\okok\\datasets"
dirresult = "E:\\R\\okok\\results"
setwd(dirresult)
getwd()
## read the datasets ##
miRNAexpression <- read.csv(paste(dirdata, "\\mirna.csv",sep=""), header = T, row.names = 1, check.names = F)
mRNAexpression  <- read.csv(paste(dirdata, "\\mrn.csv",sep=""),  header = T, row.names = 1, check.names = F)
mRNAregulator   <- read.csv(paste(dirdata, "\\result.csv",sep=""),   header = T, row.names = 1, check.names = F)
#V <- BiModule(mRNAregulator,  miRNAexpression, mRNAexpression, 0.25, 0.75)
## mainfunctions
BiModule <- function(mRNAregulator,  miRNAexpression, mRNAexpression, lambda.tol, merge.tol) {
  #record running time
  ptm <- proc.time()
  
  #index W ID to enable numerical search (more efficient than string search)
  row.id <-rownames(mRNAregulator)
  names(row.id) <-1:nrow(mRNAregulator)
  as.numeric(names(row.id))
  
  col.id <-colnames(mRNAregulator)
  names(col.id) <-1:ncol(mRNAregulator)
  as.numeric(names(col.id))
  
  dimnames(mRNAregulator)<-list(1:nrow(mRNAregulator),1:ncol(mRNAregulator))
  
  #enumerate all maximal bicliques using MICA algorithm
  Biclique <- biclique(mRNAregulator) 
  
  message("Bicluster over")
  
  #initialize the maxiaml bicliques as biclusters
  Bicluster <- bicluster(mRNAregulator, miRNAexpression)
  
  message("Bicluster over", length(Bicluster))
  
  #perform permutation test for pick out the strong bicliques
  V1 <- permutationTest(mRNAregulator, Bicluster,lambda.tol)
  
  message("permutationTest over",length(V1))
  
  #initialize co-regulatory modules
  V2 <- initmodule(mRNAregulator, V1)
  
  message("initmodule over", length(V2))
  
  #fix the co-regulatory modules assignments 
  V3 <- extendmodule(V2, mRNAregulator, merge.tol)
  
  #message("extendmodule over", length(V3))
  
  #decide the final modules
  V4 <- module(V3, mRNAregulator,  miRNAexpression, mRNAexpression, col.id, row.id) 
  
  print(sprintf("Time elapsed for BiModule: %.3f (mins)",
                (proc.time() - ptm)[3]/60))
    printModules(V4)
    return(V4) 
}
## construct weighted mRNA-regulator biclique network 
biclique <- function(W){
  #create a unweighted alternative of G, where all non-zero weights are equal to 1
  W[W>0]<-1
  
  #left, right represent the mRNAs and regulators respecitvely 
  left <- as.vector(which(W==1, arr.ind=T)[,1]) # row -> left
  right <- as.vector(which(W==1, arr.ind=T)[,2]) # col -> right
  
  edge <-data.frame(left, right)
  write.table(edge, paste(dircode, "\\bigraph.txt",sep=""), col.names = F, row.names = F, quote = F)
  write.table(edge, paste(dirresult, "\\bigraph.txt",sep=""), col.names = F, row.names = F, quote = F)
    #enumerate all maximal bicliques by employing biclique function 
  if (!require(Rcpp)) {
    install.packages('Rcpp')
  }
  library("Rcpp")
    #employing sourceCpp function to get all maximal bicliques
  sourceCpp(paste(dircode,"\\biclique.cpp",sep=""))
  biclique()
}  

## initialize the maxiaml bicliques as biclusters
bicluster <- function(W,  miRNAexpression){
    #read the biclique.txt
  bigraph <- file(paste(dirresult,"\\biclique.txt",sep=""),"r")
  
  cluster<-list(); mRNA<-c();regulator<-c()
  line <-readLines(bigraph, n=1)#read the first line of bigraph
  while(length(line) != 0 ) {
    #read mRNAs
    data1<-strsplit(line,split=" ") 
    for(i in 1:length(data1[[1]])){
      mRNA<-as.numeric(c(mRNA, data1[[1]][i]))
    }
    #read regulators
    line=readLines(bigraph,n=1)
    data2<-strsplit(line,split=" ")
    for(j in 1:length(data2[[1]])){
      regulator<-as.numeric(c(regulator,data2[[1]][j]))
    }
    line=readLines(bigraph,n=1)
    if(line=="//"){
      C <- create_v(mRNA,regulator)
      cluster <- c(cluster,list(C))
      mRNA <- c()
      regulator <-c ()
    }
    line=readLines(bigraph,n=1)
  }
  close(bigraph)
  
  #remove the the bicluster with star structures (include one regulator or one mRNA)
  cluster[sapply(cluster, function(v) (length(v$mRNA) <2 | length(v$regulator)< 2))] <- NULL
  
  return (cluster)
}

## perform a permutation test to extract those bicliques with the statistical significance ##
permutationTest<- function(W,V,lambda.tol){
  
  backgroundweights <- as.vector(W[W>0])
  
  tempV <- mclapply(V, function(v) {
    
    Observeweight <- 0; count <- 0
    
    for(i in 1:length(v$mRNA)){
      for(j in 1:length(v$regulator)){
        Observeweight <- sum(Observeweight, W[v$mRNA[i],v$regulator[j]])
      }
    }
    edgenumber <- length(v$mRNA)*length(v$regulator)
    
    #randomly sample 1000 times
    for(k in 1:1000){
      
      sampleweight <- sum(sample(backgroundweights, edgenumber, replace = F))
      
      if(Observeweight > sampleweight)
        count <- count +1
    }
    
    if(count/1000 < lambda.tol)
      v$symbol <- 1
    
    return(v)
  })
  
  tempV[sapply(tempV, function(v) (v$symbol==1))] <- NULL
  
  return (tempV)
}


## initial a module
create_v <-function(Left, Right){
  list( mRNA = Left,
        regulator = Right,
        symbol=0)
}

## Initialize TF-miRNA co-regulatory module ##
initmodule <- function(W,V){
  
  #save neighbour index for each element as a list for downstream speed-up
  nrow.list <- lapply(split(W, as.numeric(rownames(W))), 	
                      function(i) as.numeric(colnames(W)[which(i>0)]))	
  W1 <- t(W)
  ncol.list <- lapply(split(W1, as.numeric(rownames(W1))), 	
                      function(x) as.numeric(colnames(W1)[which(x>0)]))	
  
  tempV <- mclapply(V, function(v) {
    #compute D_in(G_i) and memory interactions between mRNAs and regulators at the same time
    
    countIn <- length(v$mRNA)*length(v$regulator)
    scoreIn <-0 
    for(i in 1:length(v$mRNA)){
      for(j in 1:length(v$regulator)){
        v$left <- c(v$left, v$mRNA[i])
        v$right <- c(v$right, v$regulator[j])
        scoreIn <- sum(scoreIn, W[v$mRNA[i],v$regulator[j]])
      }
    }
     scoreOutRegulator <- 0
    neighborsRegulator <- c()
    countOutRegulator <- 0
    for(i in 1:length(v$mRNA)){
      neighbors <- setdiff(nrow.list[[v$mRNA[i]]],v$regulator) #regulator neighbours for each mRNA
      neighborsRegulator <- c(neighborsRegulator, neighbors) #collect all of neighbours
    } 
    
    #
    frequenceR <- as.matrix(table(neighborsRegulator))
    frequenceR <- as.matrix(frequenceR[frequenceR[,1]>1,])
    cat("frequenceR: ", frequenceR, "\n")
    neighborsRegulator <- as.numeric(rownames(frequenceR))
    countOutRegulator <- sum(frequenceR[,1])
    if(length(neighborsRegulator)){
      for(i in 1:length(v$mRNA)){
        for(j in 1:length(neighborsRegulator)){
          scoreOutRegulator <- sum(scoreOutRegulator, W[v$mRNA[i],neighborsRegulator[j]])
        }
      }
    }
        
    
    scoreOutMRNA <-0
    neighborsMRNA <- c()
    countOutMRNA <- 0
    for(i in 1:length(v$regulator)){
      neighbors <- setdiff(ncol.list[[v$regulator[i]]],v$mRNA) #neighbour mRNAs
      neighborsMRNA <- c(neighborsMRNA, neighbors)
    }
    
    frequenceM <- as.matrix(table(neighborsMRNA))
    frequenceM <- as.matrix(frequenceM[frequenceM[,1]>1,])
    cat("frequenceM: ", frequenceM, "\n")
    neighborsMRNA <- as.numeric(rownames(frequenceM))
    countOutmRNA <- sum(frequenceM[,1])
    if(length(neighborsMRNA)){
      for(i in 1:length(v$regulator)){
        for(j in 1:length(neighborsMRNA)){
          scoreOutMRNA <- sum(scoreOutMRNA, W[neighborsMRNA[i], v$regulator[j]])
        }
      }
    }
    
    
    #compute score
    v$scoreIn <- scoreIn
    
    if((countOutRegulator+countOutMRNA) != 0)
      v$score <- scoreIn/countIn - (scoreOutRegulator + scoreOutMRNA)/(countOutRegulator + countOutMRNA)
    else 
      v$score <- scoreIn/countIn
    if(length(neighborsRegulator))
      v$neighborsRegulator <- neighborsRegulator
    else
      v$neighborsRegulator <- NULL
    if(length(neighborsMRNA))
      v$neighborsMRNA <- neighborsMRNA
    else
      v$neighborsMRNA <- NULL
    v
  })
  return (tempV)
}



## add neighbor regulator or mRNA into the biclique
add <- function(tempMRNA, tempRegulator, tempWeight, W, v, tag){
  
  scoreIn <- v$scoreIn + tempWeight 
  
  scoreOutRegulator <- 0 
  countOutRegulator <- 0
  scoreOutMRNA <- 0 
  countOutMRNA <- 0
  
  #a regulator is added  neighborsRegulator
  if(tag==1){
    
    #change the element in neighborsRegulator
    v$neighborsRegulator <- setdiff(v$neighborsRegulator, tempRegulator)
    
    if(length(v$neighborsRegulator)){
      for(i in 1:length(v$mRNA)){
        for(j in 1:length(v$neighborsRegulator)){
          scoreOutRegulator <- sum(scoreOutRegulator, W[v$mRNA[i],v$neighborsRegulator[j]])
          countOutRegulator <- countOutRegulator + 1
        }
      }
    }
    
    
    if(length(v$neighborsMRNA)){
      for(i in 1:length(v$neighborsMRNA)){
        for(j in 1:length(v$regulator)){
          scoreOutMRNA <- sum(scoreOutMRNA, W[v$neighborsMRNA[i],v$regulator[j]])
          countOutMRNA <- countOutMRNA + 1
        }
      }
    }
  }
  
  
  #a mRNA is added
  if(tag==0){ 
    #change the element in neighborsMRNA
    v$neighborsMRNA <- setdiff(v$neighborsMRNA, tempMRNA)
    
    if(length(v$neighborsMRNA)){
      for(i in 1:length(v$neighborsMRNA)){
        for(j in 1:length(v$regulator)){
          scoreOutMRNA <- sum(scoreOutMRNA, W[v$neighborsMRNA[i], v$regulator[j]])
          countOutMRNA <- countOutMRNA + 1
        }
      }
    }
    
    if(length(v$neighborsRegulator)){
      for(i in 1:length(v$mRNA)){
        for(j in 1:length(v$neighborsRegulator)){
          scoreOutRegulator <- sum(scoreOutRegulator, W[v$mRNA[i],v$neighborsRegulator[j]])
          countOutRegulator <- countOutRegulator + 1
        }
      }
    }
  }
  countIn <- length(v$left) + 1
  if((countOutRegulator+countOutMRNA) != 0)
    v$score <- scoreIn/countIn - (scoreOutRegulator + scoreOutMRNA)/(countOutRegulator + countOutMRNA)
  else 
    v$score <- scoreIn/countIn
  
  v$mRNA <- v$mRNA
  v$regulator <- v$regulator
  v$left <- c(v$left, tempMRNA)
  v$right <- c(v$right, tempRegulator)
  v$scoreIn <- scoreIn
  v$score <- v$score
  v$neighborsRegulator <- v$neighborsRegulator
  v$neighborsMRNA <- v$neighborsMRNA
  v$symbol <- 0
  return(v)
}


## merging bicliques by considering weights
extendmodule <- function(V,W,merge.tol){
  
  #decreasing order Cluster according to score 
  V1 <- V[order(sapply(V,function(v) (length(v$regulator)+length(v$mRNA))),decreasing = TRUE)]
  
  V2 <- list() #initialize new list 
  
  i <- 1 
  while(length(V1)){
    j <- 1 
    
    V2[i] <- list(V1[[j]])
    
    tempScore <- V2[[i]]$score
    
    tempMRNA      <- c()  #record mRNAs
    tempRegulator <- c()  #record regulators
    tempWeight    <- c()  #record weights between mRNAs and regulators
    
    #regulator neighbors outside of biclique 
    for(k in 1:length(V2[[i]]$mRNA)){
      neighborsRegulator <- V2[[i]]$neighborsRegulator 
      if(length(neighborsRegulator)){
        for(l in 1:length(neighborsRegulator)){
          if(W[V2[[i]]$mRNA[k],neighborsRegulator[l]]){
            tempMRNA <- c(tempMRNA, V2[[i]]$mRNA[k])
            tempRegulator <- c(tempRegulator, neighborsRegulator[l])
            tempWeight <- c(tempWeight, W[V2[[i]]$mRNA[k],neighborsRegulator[l]])
          }
        }
      }
    } 
    
    templength <- length(tempWeight) #record the number of regulator neighbors 
    
    #mRNA outside of biclique
    for(k in 1:length(V2[[i]]$regulator)){
      neighborsMRNA <- V2[[i]]$neighborsMRNA 
      if(length(neighborsMRNA)){
        for(l in 1:length(neighborsMRNA)){
          if(W[neighborsMRNA[l],V2[[i]]$regulator[k]]){
            tempMRNA <- c(tempMRNA, neighborsMRNA[l])
            tempRegulator <- c(tempRegulator, V2[[i]]$regulator[k])
            tempWeight <- c(tempWeight, W[neighborsMRNA[l], V2[[i]]$regulator[k]])
          }
        }
      }
    } 
    
    #add element into the biclique
    if(length(tempWeight)){
      #decreasing order to tempWeight
      sortWeight <- sort(tempWeight,decreasing = TRUE)
      #remember the original location
      orderWeight <- order(tempWeight,decreasing = TRUE)
      
      #adding
      for(k in 1:length(sortWeight)){
        if(orderWeight[k] <= templength){ #regulator neighbor is added
          tempV2 <- add(tempMRNA[orderWeight[k]], tempRegulator[orderWeight[k]], sortWeight[k], W, V2[[i]], 1)
          cat("Regulator tempV2$score: ", tempV2$score, " k: ", k, "\n")
          if(tempV2$score > tempScore ){
            tempScore <- tempV2$score 
            V2[i] <- list(tempV2)  #update the module
          }
          else
            break
        }
        else{ #mRNA neighbor is added
          tempV2 <-add(tempMRNA[orderWeight[k]], tempRegulator[orderWeight[k]], sortWeight[k], W, V2[[i]], 0)
          cat("mRNA tempV2$score: ", tempV2$score, " u: ", k, "\n")
          if(tempV2$score > tempScore){
            tempScore <- tempV2$score
            V2[i] <- list(tempV2)  #update the module
          }
          else
            break
        }
      }
    }
    
    #merge process
    if(length(V1)>=2){
      for(k in j+1:(length(V1)-1)){
        temp3 <- unique(unlist(V2[[i]]$left)) #mRNA
        temp4 <- unlist(V1[[k]]$regulator) #regulator
        temp5 <- unique(unlist(V2[[i]]$right)) #regulator
        temp6 <- unlist(V1[[k]]$mRNA) #mRNA
        
        if((length(intersect(temp3,temp6))+length(intersect(temp4,temp5)))
           /min((length(temp3)+length(temp5)),(length(temp4)+length(temp6))) >= merge.tol){
          V1[[k]]$symbol <- 1
        }
      }
    }
    V1[[j]] <- NULL
    V1[sapply(V1, function(v) (v$symbol==1))] <- NULL
    i <- i+1  
    
    cat("*****length(V1)***** ", length(V1), "####i#### ", i, "\n")
  }
  return (V2)
}



#calculate the density and expression correlation
module <- function(V, W,  miRNAexpression, mRNAexpression, col.id, row.id){
  
  regulatorexpression <-  miRNAexpression 
  
  tempV<- mclapply(V, function(v) {
    
    density <- 0 
    correlation <- 0
    
    for(i in 1:length(v$left)){
      correlation <- correlation + abs(as.numeric(cor(t(mRNAexpression[v$left[i],]), t(regulatorexpression[v$right
                                                                                                           [i],]))))
    }
    
    mRNA <- unique(unlist(v$left))
    regulator <- unique(unlist(v$right))
    
    density <- 2*length(v$left)/((length(regulator)+length(mRNA))*((length(regulator)+length(mRNA))-1)) 
    
    v$mRNA <- mRNA
    v$regulator <- regulator
    v$mRNA1 <- as.character(row.id[mRNA])
    v$regulator1 <- as.character(col.id[regulator])
    
    v$left <- as.character(row.id[v$left])
    v$right <- as.character(col.id[v$right])
    v$correlation <- correlation
    v$density <- density
    v$scoreIn <- NULL
    v$score <- NULL
    v$neighborsMRNA <- NULL
    v$neighborsRegulator <- NULL
    v
  })
  return (tempV)
}


## print the co-regulatory modules ##
printModules <- function(V) {
  names(V) <- paste("M",1:length(V),sep="")
  print_helper <- function(m) {
    v <- V[[m]]
    regulator <- v$regulator1
    mRNA <- v$mRNA1
    cat(sprintf("%s (density=%s,correlation=%s): \n", 
                m, v$density,v$correlation))
    cat(regulator, "\n")
    cat(mRNA, "\n")	
    cat("\\ \n")
  }	
  tmp <- lapply(names(V), print_helper)
}

## print mRNAs in each module ##
printMRNA <- function(V) {
  names(V) <- paste("M",1:length(V),sep="")
  print_helper <- function(m) {
    v <- V[[m]]
    mRNA <- v$mRNA
    cat(mRNA, "\n")	
    
  }	
  tmp <- lapply(names(V), print_helper)
}

## print mRNAs in each module ##
printRegulator <- function(V) {
  names(V) <- paste("M",1:length(V),sep="")
  print_helper <- function(m) {
    v <- V[[m]]
    regulator <- v$regulator
    cat(regulator, "\n")	
    
  }	
  tmp <- lapply(names(V), print_helper)
}


## print edges in each module ##
printEdge <- function(V) {
  names(V) <- paste("M",1:length(V),sep="")
  print_helper <- function(m) {
    v <- V[[m]]
    left <- v$left
    right <- v$right
    cat(sprintf("%s: \n", 
                m))
    cat(left, "\n")	
    cat(right,"\n")
    #cat("\\ \n")
  }	
  tmp <- lapply(names(V), print_helper)
}

regulator <- 0
mRNA <- 0
density <- c()
correlation <- c()
for(i in 1:length(V)){
  regulator <- regulator + length(V[[i]]$regulator)
  mRNA      <- mRNA + length(V[[i]]$mRNA)
  density   <- c(density, V[[i]]$density)
  correlation <- c(correlation, V[[i]]$correlation)
}

regulator <- regulator/length(V)
mRNA <- mRNA/length(V)

sum(density)/length(density)
sum(correlation)/length(correlation)

V <- BiModule(mRNAregulator,  miRNAexpression, mRNAexpression, 0.65, 0.35)

write.csv(V,"biclique-miRNA-mRNA.csv",row.names = FALSE)
########################################################################
#  step 2: MIRNA SPONGE INTERACTION NETWORK and MODULE
#########################################################################
# miRspongeR method 
##########################################################################
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

#########################################################################
# step3: function analysis
#########################################################################
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
