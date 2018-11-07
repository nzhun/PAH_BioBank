pop<-read.table("WES/PAH/UK_PAH/plink/UK.PAH.biobank.hg19.peddy.ped",comment.char = "",stringsAsFactors = F,header = 1,check.names = F,sep="\t")
## population backgroud

europeans<-c(pop$sample_id[which(pop$`ancestry-prediction`=="EUR")]);

question<-read.table("WES/PAH/UK_PAH/src/Failed.QC.ID.txt",header=F,stringsAsFactors = F)

relatedIDs<- read.table("WES/PAH/UK_PAH/src/UK.related.samples.txt",header=1,stringsAsFactors = F)
reportedNon<-read.table("WES/PAH/UK_PAH/src/report_noPAH_ID.txt",header = F,stringsAsFactors = F)
afIDs<-read.table("WES/PAH/UK_PAH/src/affected_ID.txt",header=F,stringsAsFactors = F)
afpop<-pop[which(pop$sample_id%in%afIDs$V1),]
dim(afpop)

afpop<-afpop[which(!afpop$sample_id%in%question$V1),]

dim(afpop)

afpop<-afpop[which(!afpop$sample_id%in%reportedNon$V1),]
phe<-read.csv("WES/PAH/UK_PAH/src/UK_PAH_affected_20180901.csv",header = 1,stringsAsFactors = F)
bam2ID<-read.table("WES/PAH/BAM/UK.cram2ID.txt",stringsAsFactors = F)
for(i in 1:length(phe$WGS.ID)){
  p=phe$WGS.ID[i]
  id<-grep(p,bam2ID$V2)
  if(length(id)>0){
    if(length(id)>1){print(paste(p,"error"))}
    phe$sampleID[i]<-bam2ID$V1[id]
  }else{
    print(paste(p," failed"))
  }
}
afpop<-merge(afpop,phe,by.x = "sample_id",by.y="sampleID",all.x = T)
write.table(afpop,file = "WES/PAH/UK_PAH/src/UK_PAH.affected.ped",row.names = F,sep="\t",quote=F)
dim(afpop)

# relatedIDs$GID<-unlist(lapply(1:dim(relatedIDs)[1],FUN = function(x){s=sort(relatedIDs[x,1:2]);return(paste(s[1],s[2],sep="_"))}))
# relatedIDs<-relatedIDs[!duplicated(relatedIDs$GID),]
# related_tabs<-table(c(relatedIDs$INDV1,relatedIDs$INDV2))
# rlID<-names(related_tabs[which(related_tabs>1)])
# 
# af_rlIDs<-relatedIDs[which(relatedIDs$INDV1%in% afpop[which(afpop$sample_id%in%rlID),"sample_id"] 
#                            |relatedIDs$INDV2%in% afpop[which(afpop$sample_id%in%rlID),"sample_id"]),]
# 
# mpop<-merge(afpop,af_rlIDs,by.x = "sample_id",by.y="INDV1",all.x = T)
# 
# mpop<-merge(mpop,af_rlIDs,by.x = "sample_id",by.y="INDV2",all.x = T)
# write.csv(mpop,file = "~/Desktop/af.pop.csv")