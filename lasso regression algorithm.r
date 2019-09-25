
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