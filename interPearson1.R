library(apcluster)
secSimilarity3210<-as.matrix(read.table("secSimilarity3210.txt",header=TRUE,sep=" "))
#
for(i in 1:ncol(secSimilarity3210)){
  secSimilarity3210[i,i]<-0.4
}
geneName6<-c("BRCA2","PALB2","RAD51","CHEK2","CHEK1","BARD1","BRIP1","MSH2","MSH6","MRE11A","NBN","ATM","TP53")
secSimilarity3210["BRCA1","BRCA1"]<-0.8
for(i in 1:length(geneName6)){
  secSimilarity3210[geneName6[i],geneName6[i]]<-0.6
}
#

similarity<-secSimilarity3210
m=1
totalNum<-3210
#
thresh<-0.45
threshNum<-3
increase<-0.01
factor<-50
range <- 5
minSetsize <- 2
setNum<-round(totalNum/factor)
prc <- round(range/setNum,2)*100
endSetNum <- 5
count<-0
while(setNum>=endSetNum){
  apresInterK <- apclusterK(s=similarity,K=setNum,prc=prc,maxits=1000,convits=100,lam=0.9,includeSim=TRUE, details=TRUE,nonoise=FALSE) 
  
  set<-apresInterK@clusters
  setNum<-length(apresInterK) 
  
  dirName<-paste("pearson","iter",m,"setNum",setNum,"totalNum",totalNum,"thresh",thresh,"threshNum",threshNum,sep="-")
  dir.create(dirName)
    
  newSet <- matrix()
  for(i in 1:length(set)){
    if(length(set[[i]])>=minSetsize){
      singleSet<-t(names(set[[i]]))
      singleSimilarity <- secSimilarity3210[singleSet,singleSet]
      for(j in 1:nrow(singleSimilarity)){
        if(length(which(singleSimilarity[j,]>=thresh))>=threshNum){
          newSet <- cbind(newSet,singleSet[1,j])      
        }
      }
      # 
      write.table(t(intersect(singleSet,newSet)),paste(dirName,"newSetIndividual.csv",sep="/"),
                  sep=",",append = TRUE,row.names=FALSE,col.names=FALSE,qmethod="double")
      
      rm(singleSimilarity,singleSet)
    }   
  }
  
  newSet<-newSet[1,-1]
  totalNumOld<-totalNum
  totalNum <- length(newSet)
  if(totalNumOld-totalNum<5){
    thresh<-thresh+increase
    count<-count+1
  }
  if(count==5){
    thresh<-thresh-0.03
    threshNum<-threshNum+1
    count<-0
  }
  write.table(t(newSet),paste(dirName,paste("newSetTotal",totalNum,".csv",sep="-"),sep="/"),
              sep=",",append = FALSE,row.names=FALSE,col.names=FALSE,qmethod="double")
  
  rm(similarity,apresInterK,set)
  
  similarity <-secSimilarity3210[newSet,newSet] 
  rm(newSet)
  m<-m+1
  #
  setNum<-round(totalNum/factor)
  prc <- round(range/setNum,2)*100
}







