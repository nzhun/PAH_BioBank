## UK case-control assciation
##
#setwd("~/Biocluster/")
setwd("/home/nz2274/")
source("WES/PAH/UK_PAH/script/R/UTILS.R")
#fcase="WES/PAH/UK_PAH/case_control/UK.PAH.biobank.hg19.filter.AF0.01.addregion.bed.gz"
fcase="WES/PAH/UK_PAH/case_control/UK.PAH.biobank.hg19.filter.AF0.01.vcr_xgen_union.addregion.addVQSR.bed.gz"
fctr="WES/PAH/UK_PAH/case_control/gnomad.genomes.r2.0.2.sites.NFE0.01.refGene.vcr_xgen_union.bed.MAF0.001.addregion.addDoC.corr.bed.gz"
control_dat<-read.table(fctr,comment.char = "",check.names = F,sep="\t",stringsAsFactors = F,header = 1,quote = "",strip.white = T)
case_dat<-read.table(fcase,sep="\t",stringsAsFactors = F,comment.char = "",check.names = F,header = 1,quote="")
# 
## restrict to european
#fphe<-"WES/PAH/UK_PAH/plink/UK.PAH.biobank.hg19.peddy.ped"
fphe <- "WES/PAH/UK_PAH/src/UK_PAH.affected.ped"
pop<-read.table(fphe,comment.char = "",stringsAsFactors = F,header = 1,check.names = F,sep="\t")
## population backgroud

europeans<-c(pop$sample_id[which(pop$`ancestry-prediction`=="EUR")]);
## restrict to affected

## remove related


total_case=length(europeans)
eur_dat <-case_dat[which(case_dat$proband%in%europeans),]
## filter frequency 
f=0.0001

eur_filter <- filter_case_VQSR(eur_dat,total_case,f)
eur_filter_data <-eur_filter$data
#write.csv(eur_filter$filter,"~/Desktop/case_vqsr.csv")
#eur_filter_data<-eur_filter_data[which(eur_filter_data$GQ_Ind>90),]

gnomad_filter<-filter_gnomad(control_dat,f)
gnomad_filter_data <-gnomad_filter$data
#write.csv(gnomad_filter$filter,"~/Desktop/gnomad_vqsr.csv")


cnt<-c(dim( eur_filter_data)[1]/total_case,
       dim( eur_filter_data[grep("^synony", eur_filter_data$ExonicFunc.refGene),])[1]/total_case,
       dim( eur_filter_data[grep("^nonsynony", eur_filter_data$ExonicFunc.refGene),])[1]/total_case,
       dim(eur_filter_data[grep("non",eur_filter_data$ExonicFunc.refGene,invert = T),])[1]/total_case) 

total_control<-7509
cnt2<-c(sum( gnomad_filter_data$AC_NFE)/total_control,
       sum( gnomad_filter_data[grep("^synony", gnomad_filter_data$ExonicFunc.refGene),"AC_NFE"])/total_control,
       sum( gnomad_filter_data[grep("^nonsynony", gnomad_filter_data$ExonicFunc.refGene),"AC_NFE"])/total_control,
       sum(gnomad_filter_data[grep("non",gnomad_filter_data$ExonicFunc.refGene,invert = T),"AC_NFE"])/total_control) 

write.csv(eur_filter$filter,"WES/PAH/UK_PAH/script/R/eur.dat.csv")

write.csv(gnomad_filter$filter,"WES/PAH/UK_PAH/script/R/gnomad.dat.csv")
## rest

#print(outname)
# lib_exp<-load_dataset();
# exac<-load_exac()
# sets<-load_genesets(lib_exp,exac)
genesets="ALL";
group="ALL--before"
All_genes_binom<-c()
sbinom<-test_geneset(total_case,total_control,eur_filter_data,gnomad_filter_data,genesets,1)
if(!is.null(sbinom)){
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  sbinom$Genesets<-paste(length(genesets),group)
  All_genes_binom<-rbind(All_genes_binom,sbinom)
}
write.csv(All_genes_binom,"WES/PAH/UK_PAH/case_control/result/PAH_UK_affected_hg19_VQSR.general.csv")

mdat<-merge(gnomad_filter_data,eur_filter_data)
write.table(mdat,file="WES/PAH/UK_PAH/case_control/merged_case_gnomad.bed",sep="\t",row.names=F,quote=F)

region="ALL"
process(total_case,total_control,eur_filter_data,paste("WES/PAH/UK_PAH/case_control/result/PAH-UK-hg19_affected_VQSR.",region,sep=""),gnomad_filter_data)

#x<-read.table("WES/PAH/UK_PAH/case_control/gnomad.lack.hg19.anno.addRegion.bed",sep="\t",stringsAsFactors = F,comment.char = "",check.names = F,header = 1,quote="")