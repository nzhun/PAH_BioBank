qcase<-function(data){
  return(dim(unique(data[,c("proband","Gene.refGene")]))[1])
}


test_geneset<-function(total_case,total_control,xGen_PAH,chr_vars,genesets,correct){
  keys<-c("ExonicFunc.refGene","VarClass")
  index1<-grep("ExonicFunc.refGene",names( xGen_PAH));
  if(length(index1)<1){index1<-grep("VarClass",names(xGen_PAH))};
  index2<-grep("^Func.refGene",names(xGen_PAH));
  if(length(index2)<1){index2<-grep("VarFunc",names(xGen_PAH))};
  
  ct_index1<-grep("ExonicFunc.refGene",names( chr_vars));
  if(length(ct_index1)<1){index1<-grep("VarClass",names(chr_vars))};
  ct_index2<-grep("^Func.refGene",names(chr_vars));
  if(length(ct_index2)<1){index2<-grep("VarFunc",names(chr_vars))};
  
  types<-unique(c(xGen_PAH[,index1],chr_vars[,ct_index1]))
  funs<-unique(c(xGen_PAH[,index2],chr_vars[,ct_index2]))
  
  lof_class<-types[grep("^stop|^frame",types)]  #c("stopgain","frameshiftinsertion","frameshiftdeletion","stoploss","frameshift substitution","frameshift_deletion","frameshift_insertion")
  lof_fun<-funs[grep("splic",funs)]
  binom<-c()
  rs<-c()
  subcase<-c()
  subcontrol<-c()
  #  print (gene)
  if(length(genesets)>1){
    index_case<-which(xGen_PAH$Gene.refGene%in%genesets)
    index_control<-which(chr_vars$Gene.refGene%in%genesets)
    
    if(length(index_case)<1){len_case=0}else{subcase<-xGen_PAH[index_case,] }
    if(length(index_control)<1){len_control=0;}else{  subcontrol<-chr_vars[index_control,]}
  }else{
    subcase<-xGen_PAH
    subcontrol<-chr_vars
  }
  ## synonymous
  #if(length(index_case)<1){next}
  
  if(is.null(subcase)||is.null(subcontrol)){return()}
  
  if(dim(subcase)[1]>0 && length(index1)>0){
    #   print(index1)
    #  print(names(subcase)[index1])
    len_case=qcase(subcase[grep("^synonymous",subcase[,index1]),c("proband","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0 && length(ct_index1)>0){
    len_control=sum(subcontrol$AC_NFE[grep("^synonymous",subcontrol[,ct_index1])])
  }
  test_syn<-Btest(len_case,total_case,len_control,total_control,correct) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-c("SYN",len_case,len_control,test_syn$p.value,test_syn$estimate)
  
  ### LGD
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[which(subcase[,index1] %in% lof_class |subcase[,index2]%in%lof_fun),c("proband","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0){
    len_control=sum(subcontrol$AC_NFE[which(subcontrol[,ct_index1] %in%lof_class |subcontrol[,ct_index2]%in%lof_fun)])
  }
  test_lof<-Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  rs<-rbind(rs,c("LGD",len_case,len_control,test_lof$p.value,test_lof$estimate))
  
  
  # mis
  if(dim(subcase)[1]>0){
    len_case=qcase(subcase[grep("^nonsynonymous",subcase[,index1]),c("proband","Gene.refGene")]);
  }
  if(dim(subcontrol)[1]>0){
    len_control=sum(subcontrol$AC_NFE[grep("^nonsynonymous",subcontrol[,ct_index1])])
  }
  test_mis<- Btest(len_case,total_case,len_control,total_control,correct)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
  #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
  rs<-rbind(rs,c("MIS",len_case,len_control,test_mis$p.value,test_mis$estimate))
  
  
  ### missense D (CADD>=25)
  coln<-grep("CADD13.*ph",names(subcase),ignore.case = T)[1]
  
  coln2<-grep("CADD13.*ph",names(subcontrol),ignore.case = T)[1]
  print (coln)
  print(coln2)
  for(cadd in seq(5,35,by = 5)){
    if(dim(subcase)[1]>0){
      len_case=qcase(subcase[intersect(grep("^nonsynonymous",subcase[,index1]) ,  
                                       which(as.numeric(subcase[,coln])>=cadd)),c("proband","Gene.refGene")]);
    }
    if(dim(subcontrol)[1]>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]), which(as.numeric(subcontrol[,coln2])>=cadd))])
    }
    test_cadd<- Btest(len_case,total_case,len_control,total_control,correct)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
    rs<-rbind(rs,c(paste("CADD",cadd,sep=""),len_case,len_control,test_cadd$p.value,test_cadd$estimate))
  }

  
  
  ## missense D( mcap>=0.05)
  for(thred in seq(0.025,0.2,by = 0.025)){
    if(length(subcase)>0){
      index<-intersect(grep("^nonsynonymous",subcase[,index1]), which(as.numeric(subcase$MCAP)>=thred))
      len_case<-qcase(subcase[index,]);
      #len_case=dim(unique(subcase$proband[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=0.05)),c("proband","Gene.refGene")]))[1];
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$MCAP)>=thred))])
    }
    test_mcap<- Btest(len_case,total_case,len_control,total_control,correct)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-rbind(rs,c(paste("mcap>=",thred,sep=""),len_case,len_control,test_mcap$p.value,test_mcap$estimate))
  }
 
  
  for(thred in seq(0.4,0.9,by = 0.05)){
    if(length(subcase)>0){
      
      len_case=qcase((subcase[intersect(grep("^nonsynonymous",subcase[,index1]), which(as.numeric(subcase$REVEL)>=thred)),]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol[,ct_index1]),which( as.numeric(subcontrol$REVEL)>=thred))])
    }
    #len_case=length(unique(subcase$proband[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
    #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
    test_revel<- Btest(len_case,total_case,len_control,total_control,correct) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-rbind(rs,c(paste("DREVEL>=",thred,sep=""),len_case,len_control,test_revel$p.value,test_revel$estimate))
  }
  
  binom<-rs
  ## missense D +LOF  
  #}
  #print("done")
  binom<-data.frame(as.matrix(binom),stringsAsFactors = F,check.names = F,row.names = NULL)
  
  names(binom)<-c("Type","N_case","N_control","Pvalue","OR" )
  
  return(binom)
}




process<-function(total_case,total_control,xGen_PAH,outname,chr_vars){
  #chr=8;
  print(outname)
  All_genes_binom<-c()
  types<-unique(c(xGen_PAH$ExonicFunc.refGene,chr_vars$ExonicFunc.refGene))
  
  lof_class<-types[grep("^stop|^frame",types)]  #c("stopgain","frameshiftinsertion","frameshiftdeletion","stoploss","frameshift substitution","frameshift_deletion","frameshift_insertion")
  lof_fun<-"splicing"
  binom<-c()
  nms<-c()
  
  for(gene in unique(c(chr_vars$Gene.refGene,xGen_PAH$Gene.refGene))){ #[which(xGen_PAH$`#CHROM`==chr)]
    rs<-c()
    subcase<-c()
    subcontrol<-c()
    #  print (gene)
    index_case<-which(xGen_PAH$Gene.refGene==gene)
    index_control<-which(chr_vars$Gene.refGene==gene)
    if(length(index_case)<1){len_case=0}else{subcase<-xGen_PAH[index_case,] }
    if(length(index_control)<1){len_control=0;}else{  subcontrol<-chr_vars[index_control,]}
    index_cadd_case<-grep("cadd13_phred",names(xGen_PAH),ignore.case = T)[1]
    index_cadd_ctr<-grep("cadd13_phred",names(chr_vars),ignore.case = T)[1]
    
    
    
    keys<-c("ExonicFunc.refGene","VarClass")
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[grep("^synonymous",subcase$ExonicFunc.refGene)]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[grep("^synonymous",subcontrol$ExonicFunc.refGene)])
    }
    test_syn<-Btest(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(gene,len_case,len_control,test_syn$p.value,test_syn$estimate)
    nms<-c("Gene","syN_case","syN_control","syn_Pvalue","syn_OR")
    
    ### LGD
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    }
    test_lof<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_lof$p.value,test_lof$estimate)
    nms<-c(nms,"lof_case","lof_control","lof_Pvalue","lof_OR")
    
    ### missense D (CADD>=25)
    if(length(index_cadd_case)>0 && length(index_cadd_ctr)>0){
      if(length(subcase)>0){
        len_case=length(unique(subcase$proband[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,  which(as.numeric(subcase[,index_cadd_case])>=20))]));
      }
      if(length(subcontrol)>0){
        len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene), which(as.numeric(subcontrol[,index_cadd_ctr])>=20))])
      }
      test_cadd<- Btest(len_case,total_case,len_control,total_control)  #call_pvalue(len_case,total_case,len_control,total_control) #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
      #test_cadd.estimate<-call_enrichment(len_case,len_control,total_case,total_control) 
      rs<-c(rs,len_case,len_control,test_cadd$p.value,test_cadd$estimate)
      nms<-c(nms,"cadd20_case","cadd20_control","cadd20_Pvalue","cadd20_OR")
      
      ## missense D +LOF
      if(length(subcase)>0){
        len_case=length(unique(subcase$proband[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene) ,
                                                               which(as.numeric(subcase[,index_cadd_case])>=20)), which(
                                                                 subcase$ExonicFunc.refGene %in% lof_class |
                                                                   subcase$Func.refGene%in%lof_fun)
        )])) ;
      }
      if(length(subcontrol)>0){
        
        len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) ,
                                                      which(as.numeric(subcontrol[,index_cadd_case])>=20)),which( 
                                                        subcontrol$ExonicFunc.refGene%in%lof_class |
                                                          subcontrol$Func.refGene%in%lof_fun))])
      }
      test_cadd_lof<-Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
      rs<-c(rs,len_case,len_control,test_cadd_lof$p.value,test_cadd_lof$estimate)
      nms<-c(nms,"cadd20_lof_case","cadd20_lof_control","cadd20_lof_Pvalue","cadd20_lof_OR")
      
      
      ### missense D (CADD>=25 and metasvm ==d D)
      # if(length(subcase)>0){
      #   len_case=length(unique(subcase$proband[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),which(as.numeric(subcase[,index_cadd_case])>=20  & subcase$MetaSVM_pred=="D" ))]));
      # }
      # if(length(subcontrol)>0){
      #   len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene) , which(as.numeric(subcontrol[,index_cadd_ctr])>=20 & subcontrol$MetaSVM_pred=="D"))])
      # }
      # test_cadd_meta<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
      # rs<-c(rs,len_case,len_control,test_cadd_meta$p.value,test_cadd_meta$estimate)
      # 
      ## missense D +LOF
      # if(length(subcase)>0){
      #   len_case=length(unique(subcase$proband[c(intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene),
      #                                                          which((as.numeric(subcase[,index_cadd_case])>=20  & subcase$MetaSVM_pred=="D"))) ,
      #                                                which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
      #   
      # }
      # if(length(subcontrol)>0){
      #   len_control=sum(subcontrol$AC_NFE[c(intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol[,index_cadd_ctr])>=20   & subcontrol$MetaSVM_pred=="D"))
      #                                       ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
      # }
      # 
      # test_cadd_meta_lof<- Btest(len_case,total_case,len_control,total_control)   #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
      # rs<-c(rs,len_case,len_control,test_cadd_meta_lof$p.value,test_cadd_meta_lof$estimate)
      # nms<-c(nms,"cadd20_lof_case","cadd20_lof_control","cadd20_lof_Pvalue","cadd20_lof_OR")
    }
    
    
    ## missense D( mcap>=0.05)
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$MCAP)>=0.05))]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=0.05))])
    }
    test_mcap<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_mcap$p.value,test_mcap$estimate)
    nms<-c(nms,"mcap0.05_case","mcap0.05_control","mcap0.05_Pvalue","mcap0.05_OR")
    ## missense D +LOF
    #len_case=length(unique(subcase$proband[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$MCAP)>=0.05)
    #                                                  |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun)]) );
    #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$MCAP)>=0.05)
    #                                       | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                             which((as.numeric(subcase$MCAP)>=0.05  ))) ,
                                                   which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
    }
    if(length(subcontrol)>0){
      
      len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$MCAP)>=0.05))
                                          ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
    }
    test_mcap_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_mcap_lof$p.value,test_mcap_lof$estimate)
    nms<-c(nms,"mcap0.05_lof_case","mcap0.05_lof_control","mcap0.05_lof_Pvalue","mcap0.05_lof_OR")
    
    ## missense D( revel>=0.75)
    ## missense D +LOF
    if(length(subcase)>0){
      
      len_case=length(unique(subcase$proband[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.5))]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))])
    }
    #len_case=length(unique(subcase$proband[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5)]));
    #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
    test_revel<- Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_revel$p.value,test_revel$estimate)
    nms<-c(nms,"revel0.5_case","revel0.5_control","revel0.5_Pvalue","revel0.5_OR")
    
    ## missense D +LOF
    #len_case=length(unique(subcase$proband[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.5) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
    #len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                             which((as.numeric(subcase$REVEL)>=0.5  ))) ,
                                                   which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
    }
    if(length(subcontrol)>0){
      
      len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.5))
                                          ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
    }
    test_revel_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_revel_lof$p.value,test_revel_lof$estimate)
    nms<-c(nms,"revel0.5_lof_case","revel0.5_lof_control","revel0.5_lof_Pvalue","revel0.5_lof_OR")
    
    ## missense D(mcap>=0.75)
    #len_case=length(unique(subcase$proband[which(subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75)]));
    #len_control=sum(subcontrol$AC_NFE[which(subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.5)])
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[intersect(grep("^nonsynonymous",subcase$ExonicFunc.refGene), which(as.numeric(subcase$REVEL)>=0.75))]));
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[intersect(grep("^nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.75))])
    }
    
    test_revel2<- Btest(len_case,total_case,len_control,total_control) # fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_revel2$p.value,test_revel2$estimate)
    nms<-c(nms,"revel0.75_case","revel0.75_control","revel0.75_Pvalue","revel0.75_OR")
    ## missense D +LOF
    # len_case=length(unique(subcase$proband[which((subcase$ExonicFunc.refGene=="nonsynonymous_SNV" &  as.numeric(subcase$REVEL)>=0.75) |subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun) ]));
    #  len_control=sum(subcontrol$AC_NFE[which((subcontrol$ExonicFunc.refGene=="nonsynonymous SNV" & as.numeric(subcontrol$REVEL)>=0.75) | subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun)])
    if(length(subcase)>0){
      len_case=length(unique(subcase$proband[c(intersect(grep("nonsynonymous",subcase$ExonicFunc.refGene),
                                                             which((as.numeric(subcase$REVEL)>=0.75  ))) ,
                                                   which(subcase$ExonicFunc.refGene %in% lof_class |subcase$Func.refGene%in%lof_fun))]) );
    }
    if(length(subcontrol)>0){
      len_control=sum(subcontrol$AC_NFE[c(intersect(grep("nonsynonymous",subcontrol$ExonicFunc.refGene),which( as.numeric(subcontrol$REVEL)>=0.75))
                                          ,which(subcontrol$ExonicFunc.refGene%in%lof_class |subcontrol$Func.refGene%in%lof_fun))])
    }
    test_revel2_lof<- Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
    rs<-c(rs,len_case,len_control,test_revel2_lof$p.value,test_revel2_lof$estimate)
    nms<-c(nms,"revel0.75_lof_case","revel0.75_lof_control","revel0.75_lof_Pvalue","revel0.75_lof_OR")
    
    binom<-rbind(binom,rs)
    ## missense D +LOF  
  }
  #print("done")
  binom<-data.frame(as.matrix(binom),stringsAsFactors = F,check.names = F,row.names = NULL)
  print(dim(binom))
  names(binom)<-nms
  
  # names(binom)<-c("Gene","syN_case","syN_control","syn_Pvalue","syn_OR",
  #                 "LOF_case","LOF_control","LOF_Pvalue","LOF_OR",
  #                 "CaddN_case","CaddN_control","Cadd_Pvalue","Cadd_OR",
  #                 "Cadd_LOF_N_case","Cadd_LOF_N_control","Cadd_LOF_Pvalue","Cadd_LOF_OR",
  #                 "CaddN_Meta_case","CaddN_Meta_control","Cadd_Meta_Pvalue","Cadd_Meta_OR",
  #                 "Cadd_Meta_LOF_N_case","Cadd_Meta_LOF_N_control","Cadd_Meta_LOF_Pvalue","Cadd_Meta_LOF_OR",
  #                 "MCAP_N_case","MCAP_N_control","MCAP_Pvalue","MCAP_OR",
  #                 "MCAP_LOF_N_case","MCAP_LOF_N_control","MCAP_LOF_Pvalue","MCAP_LOF_OR",
  #                 "REVEL_N_case","REVEL_N_control","REVEL_Pvalue","REVEL_OR",
  #                 "REVEL_LOF_N_case","REVEL_LOF_N_control","REVEL_LOF_Pvalue","REVEL_LOF_OR",
  #                 "REVEL2_N_case","REVEL2_N_control","REVEL2_Pvalue","REVEL2_OR",
  #                 "REVE2L_LOF_N_case","REVEL2_LOF_N_control","REVEL2_LOF_Pvalue","REVEL2_LOF_OR"
  # )
  #    write.csv(binom,file = paste("PAH/PAH-CHD/",outname,".gnomeAD.",chr,".binom.csv",sep=""),row.names = F)
  All_genes_binom<-rbind(All_genes_binom,binom)
  
  
  
  All_genes_binom$LOF_pvalue_adjust<-p.adjust(All_genes_binom$LOF_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
  
  All_genes_binom$REVEL_LOF_Pvalue_adjust<-p.adjust(All_genes_binom$REVEL_LOF_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
  
  All_genes_binom$REVEL_Pvalue_adjust<-p.adjust(All_genes_binom$REVEL_Pvalue,method = p.adjust.methods[5],n = dim(All_genes_binom)[1])
  
  outf= paste(outname,".gnomeAD.ALL.binom.csv",sep="")
  write.csv(All_genes_binom,file =outf,row.names = F)
  
  fout=paste(outname,".gnomeAD.11PAH.binom.csv",sep="")  #"~/PAH/PAH-CHD/PAH-CHD.gnomeAD.253CHD.binom.csv"
  fgene="WES/PAH/Result/Data/source/PAH_associated11-13.txt" #"~/PAH/Result/Data/source/CHD.253GeneList.csv"
  
  genes<-read.csv(fgene,header=1,stringsAsFactors=F,check.names=F,comment.char="",quote="")
  sub<-All_genes_binom[which(All_genes_binom$Gene%in%genes$Gene),]
  allr<-c("Gene")
  if(dim(sub)[1]>0){
    for(i in seq(2,(dim(sub)[2]-3),by = 4)){
      len_case<-as.integer(sum(as.numeric(sub[,i])))
      len_control<-as.integer(sum(as.numeric(sub[,(i+1)])))
      allr<-c(allr,c(len_case,len_control))
      test<- Btest(len_case,total_case,len_control,total_control)   #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
      allr<-c(allr,c(as.numeric(test$p.value),as.numeric(test$estimate)))
    }
    allr<-c(allr,rep("-",dim(sub)[2]-length(allr)))
    #  allr<-data.frame(allr,stringsAsFactors = F)
    #names(allr)<-names(sub)
    write.csv(rbind(as.matrix(sub),allr),fout,row.names=F)
  }
  
  fout=paste(outname,".gnomeAD.253CHD.binom.csv",sep="")  #"~/PAH/PAH-CHD/PAH-CHD.gnomeAD.253CHD.binom.csv"
  fgene="WES/PAH/Result/Data/source/CHD.253GeneList.csv"
  
  genes<-read.csv(fgene,header=1,stringsAsFactors=F,check.names=F,comment.char="",quote="")
  sub<-All_genes_binom[which(All_genes_binom$Gene%in%genes$Gene),]
  allr<-c("Gene")
  if(dim(sub)[1]>0){
    for(i in seq(2,(dim(sub)[2]-3),by = 4)){
      # print(i)
      len_case<-as.integer(sum(as.numeric(sub[,i])))
      len_control<-as.integer(sum(as.numeric(sub[,(i+1)])))
      allr<-c(allr,c(len_case,len_control))
      test<-Btest(len_case,total_case,len_control,total_control)  #fisher.test(matrix(c(len_case,total_case-len_case,len_control,total_control-len_control),nrow = 2,byrow = T))
      allr<-c(allr,c(as.numeric(test$p.value),as.numeric(test$estimate)))
    }
    allr<-c(allr,rep("-",dim(sub)[2]-length(allr)))
    #  allr<-data.frame(allr,stringsAsFactors = F)
    #  names(allr)<-names(sub)
    write.csv(rbind(as.matrix(sub),allr),fout,row.names=F)   
    print(fout)
    #  }
    
  }
  
  
  
  All_genes_binom<-read.csv(outf,header = 1,check.names = T,stringsAsFactors = F)
  
  
  
  
  pdf(paste(outname,".gnomeAD.QQ.pdf",sep=""),width = 5,height = 5)
  par(mai=c(0.5,0.5,0,0),mar=c(4,4,2,1))
  if(!is.null(All_genes_binom) && dim(All_genes_binom)[1]>0){
    x<--log10((dim(All_genes_binom)[1]:1)/dim(All_genes_binom)[1])
    y<--log10(sort(as.numeric(All_genes_binom$syn_Pvalue),decreasing = T))
    if(length(y)>0){
      plot(x,y,main="SYN QQ plot",xlab="-log10 (expected p-value)",ylab="-log10 (observed p-value)",xlim=c(0,max(y,x)),ylim=c(0,max(y,x)))
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
    }
    x<--log10((dim(All_genes_binom)[1]:1)/dim(All_genes_binom)[1])
    y<--log10(sort(as.numeric(All_genes_binom$LOF_Pvalue),decreasing = T))
    if(length(y)>0){
      
      plot(x,y,main="LGD QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",pch=20)
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      #  points( max(x) ,max(y),pch=5,col="red")
      
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$LOF_Pvalue==min(as.numeric(All_genes_binom$LOF_Pvalue))],pos = 2,col="red")
    }
    
    
    x<--log10((dim(All_genes_binom)[1]:1)/dim(All_genes_binom)[1])
    # y<--log10(sort(as.numeric(All_genes_binom$LOF_pvalue_adjust),decreasing = T))
    # plot(x,y,main="LGD QQ plot",xlab="-log10(expected)",ylab="-log10(observed_adjusted p-value)")
    # abline(0,1,lty=2,col="gray")
    # abline(h=-log10(0.05),lty=2,col="gray")
    # 
    # text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$LOF_pvalue_adjust==min(as.numeric(All_genes_binom$LOF_pvalue_adjust))],pos = 2,col="red")
    
    
    
    y<--log10(sort(as.numeric(All_genes_binom$Cadd_Pvalue),decreasing = T))
    if(length(y)>0){
      plot(x,y,main="CADD>=20+ LGD QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)",xlim=c(0,max(y,x)),ylim=c(0,max(y,x)))
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$Cadd_Pvalue==min(as.numeric(All_genes_binom$Cadd_Pvalue))])
    }
    
    y<--log10(sort(as.numeric(All_genes_binom$Cadd_LOF_Pvalue),decreasing = T))
    if(length(y)>0){
      plot(x,y,main="CADD>=20 QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)")
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$Cadd_LOF_Pvalue==min(as.numeric(All_genes_binom$Cadd_LOF_Pvalue))],pos = 2,col="red")
    }
    
    y<--log10(sort(as.numeric(All_genes_binom$REVEL_Pvalue),decreasing = T))
    if(length(y)>0){
      plot(x,y,main="REVEL>=0.5 QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)")
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$REVEL_Pvalue==min(as.numeric(All_genes_binom$REVEL_Pvalue))],pos = 2,col="red")
    }
    
    y<--log10(sort(as.numeric(All_genes_binom$REVEL_LOF_Pvalue),decreasing = T))
    if(length(y)>0){
      plot(x,y,main="",xlab="",
           ylab="",pch=20,xlim=c(0,max(y,x)),
           ylim=c(0,max(y,x)),cex.axis=1.5,cex.lab=1.5)
      abline(0,1,lty=2,col="gray")
      abline(h=-log10(0.05/20000),lty=2,col="gray")
      title(main="REVEL>0.5 and LGD QQ plot")
      title(xlab="-log10 (expected p-value)",ylab="-log10 (observed p-value)",line = 2.5,cex.lab=1.5)
      text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$REVEL_LOF_Pvalue==min(as.numeric(All_genes_binom$REVEL_LOF_Pvalue))],
           font=3,pos = 4,col="red",cex=1.5)
      points( max(x) ,max(y),pch=19,col="red",cex=2)
    }
    # y<--log10(sort(as.numeric(All_genes_binom$REVEL_LOF_Pvalue_adjust),decreasing = T))
    # plot(x,y,main="REVEL>=0.5 +LGD QQ plot",xlab="-log10(Expected p-value)",ylab="-log10(Observed p-value)")
    # abline(0,1,lty=2,col="gray")
    # abline(h=-log10(0.05),lty=2,col="gray")
    # text(max(x),max(y),All_genes_binom$Gene[All_genes_binom$REVEL_LOF_Pvalue_adjust==min(as.numeric(All_genes_binom$REVEL_LOF_Pvalue_adjust))],pos = 2,col="red")
    
    dev.off()
  }
  
}



format_AAchange<-function(dat){
  h<-grep("AAChange",names(dat))
  if(length(h)>0){
    dat$GeneName<-unlist(lapply(dat[,h],FUN = function(x)unlist(strsplit(x,split = ":"))[1] ))
    dat$Transcript<-unlist(lapply(dat[,h],FUN = function(x)unlist(strsplit(x,split = ":"))[2] ))
    dat$Exon<-unlist(lapply(dat[,h],FUN = function(x)unlist(strsplit(x,split = ":"))[3] ))
    
    dat$NucleiotideChange<-unlist(lapply(dat[,h],FUN = function(x)unlist(strsplit(x,split = ":"))[4] ))
    dat$ProteinChange<-unlist(lapply(dat[,h],FUN = function(x)unlist(strsplit(x,split = ":"))[5] ))
    
    dat$ProteinChange<-unlist(lapply(dat[,h],FUN = function(x)unlist(strsplit(x,split = ":"))[5] ))
    
    
    h<-which(names(dat)=="AAChange.refGene")
    v<-grep("VarClass|ExonicFunc.refGene",names(dat))
    f<-grep("VarFunc|Func.refGene",names(dat))
    if(length(h)>0){
      dat$Mutation_Type<-""
      dat$Mutation_Type[grep("stop",dat[,v])]<-"LGD"
      dat$Mutation_Type[grep("^frame",dat[,v])]<-"LGD"
      dat$Mutation_Type[grep("nonframe",dat[,v])]<-"In_Frame"
      dat$Mutation_Type[grep("nonsynonymous",dat[,v])]<-"Mis"
      dat$Mutation_Type[grep("^synony",dat[,v])]<-"SYN"
      dat$Mutation_Type[which(dat[,f]=="splicing")]<-"LGD"
    } 
  }
  
  return(dat)
  
}


### number of observation in cases, number of observation in control 
Btest<-function(case,case_total,control,control_total,correct,side){
  if(missing(correct)){correct=1}
  if(missing(side)){side="two.sided"}
  case<-round(case/correct)
  if(case+control>0){
    btest<- binom.test(case,case+control,p=case_total/(control_total+case_total),alternative = side,conf.level = 0.95)
    or<-(case/control) * (control_total/case_total)
    btest<-data.frame(matrix(c(btest$p.value,or),nrow=1),check.names = F,stringsAsFactors = F);
  }else{
    btest<-data.frame(matrix(c(1,1),nrow=1),check.names = F,stringsAsFactors = F);
    
  }
  names(btest)<-c("p.value","estimate")
  return (btest)
}
call_pvalue<-function(case,case_total,control,control_total,side){
  if(missing(side)){side="two.sided"}
  #return(binom.test(case,case_total,p=control/control_total,alternative = side,conf.level = 0.95))
  return(binom.test(case,case+control,p=case_total/(control_total+case_total),alternative = side,conf.level = 0.95))
}
call_enrichment<-function(case,control,len_case,len_control,correction){
  return ((length(case)/correction)*len_control/(len_case*length(control)))
  
}

call_enrichment_n<-function(lcase,lcontrol,len_case,len_control){
  return (lcase*len_control/(len_case*lcontrol))
  
}


get_lgd<-function(dat){
  index1<-grep("stop",dat$VarClass)
  index2<-grep("^frame",dat$VarClass)
  index3<-grep("splicing",dat$VarFunc)
  return(unique(c(index1,index2,index3)))
}
get_indel<-function(dat){
  len_ref<-unlist(lapply(dat$REF,FUN = function(x) nchar(x)))
  len_alt<-unlist(lapply(dat$ALT,FUN=function(x)nchar(x)))
  if(length(len_ref)!=length(len_alt)){stop("Error: the number of ref and alt are not equal!"); }
  return(which(len_ref!=len_alt))
}
get_snv<-function(dat){
  len_ref<-unlist(lapply(dat$REF,FUN = function(x) nchar(x)))
  len_alt<-unlist(lapply(dat$ALT,FUN=function(x)nchar(x)))
  if(length(len_ref)!=length(len_alt)){stop("Error: the number of ref and alt are not equal!"); }
  return(which(len_ref==len_alt))
}



formatFreq<-function(data){
  if(length(grep("1KGfreq",names(data)))>0){
    data$`1KG.afr.freq`[which(is.na(as.numeric(data$`1KG.afr.freq`)))]=0
    data$`1KG.eas.freq`[which(is.na(as.numeric(data$`1KG.eas.freq`)))]=0
    data$`1KG.amr.freq`[which(is.na(as.numeric(data$`1KG.amr.freq`)))]=0
    data$`1KG.eur.freq`[which(is.na(as.numeric(data$`1KG.eur.freq`)))]=0
    data$`1KG.sas.freq`[which(is.na(as.numeric(data$`1KG.sas.freq`)))]=0
    data$`1KGfreq`[which(is.na(as.numeric(data$`1KGfreq`)))]=0
  }
  if(length(grep("ESPfreq",names(data)))>0){
    data$ESPfreq[which(is.na(as.numeric(data$ESPfreq)))]=0
    data$ESP.aa.freq[which(is.na(as.numeric(data$ESP.aa.freq)))]=0
    data$ESP.ea.freq[which(is.na(as.numeric(data$ESP.ea.freq)))]=0
  }
  if(length(grep("ExACfreq",names(data)))>0){
    data$ExACfreq[which(is.na(as.numeric(data$ExACfreq)))]=0
    data$ExAC.eas.freq[which(is.na(as.numeric(data$ExAC.eas.freq)))]=0
    data$ExAC.afr.freq[which(is.na(as.numeric(data$ExAC.afr.freq)))]=0
    data$ExAC.amr.freq[which(is.na(as.numeric(data$ExAC.amr.freq)))]=0
    data$ExAC.nfe.freq[which(is.na(as.numeric(data$ExAC.nfe.freq)))]=0
    data$ExAC.sas.freq[which(is.na(as.numeric(data$ExAC.sas.freq)))]=0
    data$ExAC.fin.freq[which(is.na(as.numeric(data$ExAC.fin.freq)))]=0
    data$ExAC.oth.freq[which(is.na(as.numeric(data$ExAC.oth.freq)))]=0
  }
  if(length(grep("gnomAD_Exome_AF",names(data)))>0){
    data$gnomAD_Exome_AF[which(is.na(as.numeric(data$gnomAD_Exome_AF)))]=0
    data$gnomAD_Exome_AF_ASJ[which(is.na(as.numeric(data$gnomAD_Exome_AF_ASJ)))]=0
    data$gnomAD_Exome_AF_AMR[which(is.na(as.numeric(data$gnomAD_Exome_AF_AMR)))]=0
    data$gnomAD_Exome_AF_AFR[which(is.na(as.numeric(data$gnomAD_Exome_AF_AFR)))]=0
    data$gnomAD_Exome_AF_EAS[which(is.na(as.numeric(data$gnomAD_Exome_AF_EAS)))]=0
    data$gnomAD_Exome_AF_SAS[which(is.na(as.numeric(data$gnomAD_Exome_AF_SAS)))]=0
    data$gnomAD_Exome_AF_OTH[which(is.na(as.numeric(data$gnomAD_Exome_AF_OTH)))]=0
    data$gnomAD_Exome_AF_NFE[which(is.na(as.numeric(data$gnomAD_Exome_AF_NFE)))]=0
    data$gnomAD_Exome_AF_FIN[which(is.na(as.numeric(data$gnomAD_Exome_AF_FIN)))]=0
  }
  if(length(grep("gnomAD_Genome_AF",names(data)))>0){
    data$gnomAD_Genome_AF[which(is.na(as.numeric(data$gnomAD_Genome_AF)))]=0
    data$gnomAD_Genome_AF_AFR[which(is.na(as.numeric(data$gnomAD_Genome_AF_AFR)))]=0
    data$gnomAD_Genome_AF_ASJ[which(is.na(as.numeric(data$gnomAD_Genome_AF_ASJ)))]=0
    data$gnomAD_Genome_AF_AMR[which(is.na(as.numeric(data$gnomAD_Genome_AF_AMR)))]=0
    data$gnomAD_Genome_AF_AFR[which(is.na(as.numeric(data$gnomAD_Genome_AF_AFR)))]=0
    data$gnomAD_Genome_AF_EAS[which(is.na(as.numeric(data$gnomAD_Genome_AF_EAS)))]=0
    data$gnomAD_Genome_AF_OTH[which(is.na(as.numeric(data$gnomAD_Genome_AF_OTH)))]=0
    data$gnomAD_Genome_AF_NFE[which(is.na(as.numeric(data$gnomAD_Genome_AF_NFE)))]=0
    data$gnomAD_Genome_AF_FIN[which(is.na(as.numeric(data$gnomAD_Genome_AF_FIN)))]=0
    data$gnomAD_Genome_AF_POPMAX[which(is.na(as.numeric(data$gnomAD_Genome_AF_POPMAX)))]=0
  }
  return(data)
}


formatFreq_new<-function(data){
  if(length(grep("1000g2015aug_all",names(data)))>0){
    data$`1000g2015aug_afr`[which(is.na(as.numeric(data$`1000g2015aug_afr`)))]=0
    data$`1000g2015aug_eas`[which(is.na(as.numeric(data$`1000g2015aug_eas`)))]=0
    data$`1000g2015aug_amr`[which(is.na(as.numeric(data$`1000g2015aug_amr`)))]=0
    data$`1000g2015aug_eur`[which(is.na(as.numeric(data$`1000g2015aug_eur`)))]=0
    data$`1000g2015aug_sas`[which(is.na(as.numeric(data$`1000g2015aug_sas`)))]=0
    data$`1000g2015aug_all`[which(is.na(as.numeric(data$`1000g2015aug_all`)))]=0
  }
  if(length(grep("ExAC_ALL",names(data)))>0){
    data$ExAC_ALL[which(is.na(as.numeric(data$ExAC_ALL)))]=0
    data$ExAC_EAS[which(is.na(as.numeric(data$ExAC_EAS)))]=0
    data$ExAC_AFR[which(is.na(as.numeric(data$ExAC_AFR)))]=0
    data$ExAC_AMR[which(is.na(as.numeric(data$ExAC_AMR)))]=0
    data$ExAC_NFE[which(is.na(as.numeric(data$ExAC_NFE)))]=0
    data$ExAC_SAS[which(is.na(as.numeric(data$ExAC_SAS)))]=0
    data$ExAC_FIN[which(is.na(as.numeric(data$ExAC_FIN)))]=0
    data$ExAC_OTH[which(is.na(as.numeric(data$ExAC_OTH)))]=0
  }
  if(length(grep("esp6500siv2_all",names(data)))>0){
    data$esp6500siv2_all[which(is.na(as.numeric(data$esp6500siv2_all)))]=0
    data$esp6500siv2_aa[which(is.na(as.numeric(data$esp6500siv2_aa)))]=0
    data$esp6500siv2_ea[which(is.na(as.numeric(data$esp6500siv2_ea)))]=0
  }
  if(length(grep("gnomAD_exome_ALL",names(data)))>0){
    data$gnomAD_exome_ALL[which(is.na(as.numeric(data$gnomAD_exome_ALL)))]=0
    data$gnomAD_exome_ASJ[which(is.na(as.numeric(data$gnomAD_exome_ASJ)))]=0
    data$gnomAD_exome_AMR[which(is.na(as.numeric(data$gnomAD_exome_AMR)))]=0
    data$gnomAD_exome_AFR[which(is.na(as.numeric(data$gnomAD_exome_AFR)))]=0
    data$gnomAD_exome_EAS[which(is.na(as.numeric(data$gnomAD_exome_EAS)))]=0
    data$gnomAD_exome_FIN[which(is.na(as.numeric(data$gnomAD_exome_FIN)))]=0
    data$gnomAD_exome_SAS[which(is.na(as.numeric(data$gnomAD_exome_SAS)))]=0
    data$gnomAD_exome_OTH[which(is.na(as.numeric(data$gnomAD_exome_OTH)))]=0
    data$gnomAD_exome_NFE[which(is.na(as.numeric(data$gnomAD_exome_NFE)))]=0
  }
  if(length(grep("gnomAD_genome_ALL",names(data)))>0){
    data$gnomAD_genome_ALL[which(is.na(as.numeric(data$gnomAD_genome_ALL)))]=0
    data$gnomAD_genome_AFR[which(is.na(as.numeric(data$gnomAD_genome_AFR)))]=0
    data$gnomAD_genome_ASJ[which(is.na(as.numeric(data$gnomAD_genome_ASJ)))]=0
    data$gnomAD_genome_AMR[which(is.na(as.numeric(data$gnomAD_genome_AMR)))]=0
    data$gnomAD_genome_FIN[which(is.na(as.numeric(data$gnomAD_genome_FIN)))]=0
    data$gnomAD_genome_AFR[which(is.na(as.numeric(data$gnomAD_genome_AFR)))]=0
    data$gnomAD_genome_EAS[which(is.na(as.numeric(data$gnomAD_genome_EAS)))]=0
    data$gnomAD_genome_OTH[which(is.na(as.numeric(data$gnomAD_genome_OTH)))]=0
    data$gnomAD_genome_NFE[which(is.na(as.numeric(data$gnomAD_genome_NFE)))]=0
    # data$gnomAD_genome_POPMAX[which(is.na(as.numeric(data$gnomAD_genome_POPMAX)))]=0
  }
  return(data)
}



filter_allfreq <- function(data,freq_avg,freq_max){
  #  freq <- 0.001
  #  freq2 <- 0.001
  l<-length(grep("ExACfreq",names(data)))
  
  if(l<1) return( data);
  
  data <- data[which(na.pass(as.numeric(data$ExACfreq)< freq_avg)
                     &na.pass(as.numeric(data$ExAC.amr.freq)< freq_max)
                     &as.numeric(data$ExAC.afr.freq)< freq_max
                     &as.numeric(data$ExAC.nfe.freq)< freq_max
                     &as.numeric(data$ExAC.sas.freq)< freq_max
                     &as.numeric(data$ExAC.eas.freq)< freq_max
                     &as.numeric(data$ExAC.oth.freq)< freq_max
                     #  & as.numeric(data$`1KGfreq`) < freq2
                     #  &as.numeric(data$ESPfreq)< freq2
                     #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  #  data <- data[which( as.numeric(data$AC)< 25
  #                     &as.numeric(data$AB)>0.2
  # ),]
  if(length(grep("Mappability",names(data)))>0){ 
    data <- data[which(data$Mappability==1),]
  }
  return(data)
}


load_lof_CADDscore<-function(gene,dataset){
  ### /home/local/ARCS/nz2274/Application/psap/lookups/full.lof.pCADD.gencodeV19.allsites.txt.gz
  record<-0
  len<-length(which(dataset$V1==gene))
  if(len>0)  record<-dataset[which(dataset$V1==gene),3]
  return(record)  
}

get_popscore<- function(gene,score,dataset){
  ##/home/local/ARCS/nz2274/Application/psap/lookups/full.het.pCADD.gencodeV19.allsites.txt.gz
  popscore<-1;
  scale<-seq(0,70,0.05)
  len<-length(which(dataset$V1==gene))
  if(len>0){
    vec<-as.numeric(dataset[which(dataset$V1==gene),2:dim(dataset)[2]])
    popscore<-vec[findInterval(score,scale)+1]
  }
  return(popscore)
}

filter_allfreq_new <- function(data,freq_avg,freq_max){
  
  data <- data[which(na.pass(as.numeric(data$ExAC_ALL)< freq_avg)
                     &na.pass(as.numeric(data$ExAC_AMR)< freq_max)
                     &as.numeric(data$ExAC_AFR)< freq_max
                     &as.numeric(data$ExAC_NFE)< freq_max
                     &as.numeric(data$ExAC_FIN)< freq_max
                     &as.numeric(data$ExAC_SAS)< freq_max
                     &as.numeric(data$ExAC_EAS)< freq_max
                     &as.numeric(data$ExAC_OTH)< freq_max
                     &as.numeric(data$gnomAD_genome_EAS)<freq_max
                     &as.numeric(data$gnomAD_genome_NFE)<freq_max
                     &as.numeric(data$gnomAD_genome_FIN)<freq_max
                     &as.numeric(data$gnomAD_genome_OTH)<freq_max
                     &as.numeric(data$gnomAD_genome_ASJ)<freq_max
                     &as.numeric(data$gnomAD_genome_AMR)<freq_max
                     &as.numeric(data$gnomAD_genome_ALL)<freq_avg
                     &as.numeric(data$gnomAD_exome_ALL)<freq_avg
                     &as.numeric(data$gnomAD_exome_EAS)<freq_max
                     &as.numeric(data$gnomAD_exome_NFE)<freq_max
                     &as.numeric(data$gnomAD_exome_FIN)<freq_max
                     &as.numeric(data$gnomAD_exome_OTH)<freq_max
                     &as.numeric(data$gnomAD_exome_ASJ)<freq_max
                     &as.numeric(data$gnomAD_exome_AMR)<freq_max
                     #  & as.numeric(data$`1KGfreq`) < freq2
                     #  &as.numeric(data$ESPfreq)< freq2
                     #  &as.numeric(data$gnomAD_Genome_AF)< freq2
  ),]
  # data <- data[which( as.numeric(data$AC)< 25
  #                      &as.numeric(data$AB)>0.2
  #  ),]
  # if(length(grep("Mappability",names(data)))>0){ 
  #   indexs<-which(data$Mappability==1)
  #   if(length(indexs)>0){
  #     data <- data[indexs,]
  #   }
  # }
  # if(length(grep("genomicSuperDups",names(data)))>0){
  #   index<-grep("Score",data$genomicSuperDups)
  #   as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))
  #   dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-"))[2])))>0.9)]
  #   if(length(dup_indexs)>0){
  #     data <- data[-dup_indexs,]
  #   }
  # }
  return(data)
}



load_dataset<-function(){
  resources="Resources/"
  # fexac<-paste(resources,"ExAC_r03_0316_z_pli_rec_null_data.bed",sep="")
  fheart<-paste(resources,"GeneExpression/mousebrain_Heart.csv",sep="")
  flung<-paste(resources,"GeneExpression/Lung_rank_RNAseq Asselin-Labat-GPL13112_human.csv",sep="")
  fdiaphram<-paste(resources,"GeneExpression/diaphragm_rank.csv",sep="")
  # exac<-read.table(fexac,stringsAsFactors = F,header = 1,check.names = F,comment.char = "")
  heart_exp<-read.csv(fheart,sep = ",",header = 1,check.names = F)
  lung_exp<-read.csv(flung,sep=",",header = 1,check.names = F)
  diaphragm_exp<-read.csv(fdiaphram,sep=",",header = 1,check.names = F)
  
  
  return(list(heart_exp=heart_exp,lung_exp=lung_exp,diaph_exp=diaphragm_exp))
}
load_exac<-function(){
  resources="Resources/"
  fexac<-paste(resources,"ExAC_r03_0316_z_pli_rec_null_data.bed.gz",sep="")
  exac<-read.table(fexac,stringsAsFactors = F,header = 1,check.names = F,comment.char = "")
  return (exac)
}
load_roadmap_RKPM<-function(){
  resources="~/server/Resources/"
  
  file<-paste(resources,"GeneExpression/57epigenomes.RPKM.pc.csv.gene.txt",sep="")
  exp<-read.table(file,header = 1,sep = "\t",check.names = F,stringsAsFactors = F,row.names = NULL)
  #flung<-paste(resources,"GeneExpression/Lung_rank_RNAseq Asselin-Labat-GPL13112_human.csv",sep="")
  fdiaphram<-paste(resources,"GeneExpression/diaphragm_rank.csv",sep="")
  aorta<-"E065"
  lung<-"E096"
  heart_exp<-exp[,c("GeneName",aorta)]
  lung_exp<-exp[,c("GeneName",lung)]
  diaphragm_exp<-read.csv(fdiaphram,sep=",",header = 1,check.names = F,strip.white = T)
  
  return(list(heart_exp=heart_exp,lung_exp=lung_exp,diaph_exp=diaphragm_exp))
}



format_consensus<-function(data){
  nms<-names(data)
  nms[which(nms=="Gene.refGene")]="GeneName"
  nms[which(nms=="ExonicFunc.refGene")]="VarClass"
  nms[which(nms=="Func.refGene")]="VarFunc"
  nms[which(nms=="ExAC_ALL")]="ExACfreq"
  nms[which(nms=="ExAC_AMR")]="ExAC.amr.freq"
  nms[which(nms=="ExAC_EAS")]="ExAC.eas.freq"
  nms[which(nms=="ExAC_NFE")]="ExAC.nfe.freq"
  nms[which(nms=="ExAC_AFR")]="ExAC.afr.freq"
  nms[which(nms=="ExAC_FIN")]="ExAC.fin.freq"
  nms[which(nms=="1000g2015aug_afr")]="1KG.afr.freq"
  nms[which(nms=="1000g2015aug_all")]="1KG.all.freq"
  nms[which(nms=="1000g2015aug_amr")]="1KG.amr.freq"
  nms[which(nms=="1000g2015aug_eas")]="1KG.eas.freq"
  nms[which(nms=="1000g2015aug_eur")]="1KG.eur.freq"
  nms[which(nms=="1000g2015aug_sas")]="1KG.sas.freq"
  nms[which(nms=="AAChange.refGene")]="AAChange"
  #  nms[which(nms=="CADD13_phred")]="CADDphred"
  nms[which(nms=="CADD_phred")]="CADDphred"
  nms[which(nms=="CADDraw")]="CADD13_raw"
  # nms[which(nms=="VarFunc")]="Func.refGene"
  nms[which(nms=="MetaLR_pred")]="MetaLRprd"
  nms[which(nms=="MetaSVM_pred")]="MetaSVMprd"
  nms[which(nms=="MetaSVMscr")]="MetaSVM_score"
  nms[which(nms=="MutAprd")]="MutationAssessor_pred"
  nms[which(nms=="MutAscr")]="MutationAssessor_score"
  nms[which(nms=="MutTscr")]="MutationTaster_score" 
  nms[which(nms=="MutTprd")]="MutationTaster_pred" 
  nms[which(nms=="Polyphen2_HVAR_pred")]="PP2.hvar.prd" 
  nms[which(nms=="Polyphen2_HDIV_pred")]="PP2.hdiv.prd" 
  #nms[which(nms=="proband")]="proband_GT"
  if(length(which(nms=="proband_GT"))==0){
    data$proband_GT<-data$proband
  }
  data$proband<-data$proband
  #nms[which(nms=="proband")]="proband"
  
  names(data)<-nms
  return(data)
  
}

gene_exp_diaphragm<-function(data,lib){
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")
    index<-grep(data[id,gid],lib$gene)[1]
    if(length(index)>0){
      data$diaphragm_rank[id]=lib$rank[index]
      data$diaphragm_E115[id]=lib$E11.5[index]
    }else{
      data$diaphragm_rank[id]="."
      data$diaphragm_E115[id]="."
      
    }
  }
  return(data)
}

gene_exp_Heart<-function(data,lib){
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")
    index<-grep(data[id,gid],lib$human.External.Gene.Name)[1]
    if(length(index)>0){
      data$HEART_EXP[id]=lib$e14.5_rank[index]
      #  data$diaphragm_E115[id]=lib$E11.5[index]
    }else{
      data$HEART_EXP[id]="."
      #data$diaphragm_E115[id]="."
      
    }
  }
  print("Heart Expression, the higher the more related")
  return(data)
}



gene_exp_Lung<-function(data,lib){ ## reversed
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")[1]
    
    index<-grep(data[id,gid],lib$human_gene)[1]
    if(length(index)>0){
      data$LUNG_EXP[id]=lib$`Control-Stroma rank`[index]
      #  data$diaphragm_E115[id]=lib$E11.5[index]
    }else{
      data$LUNG_EXP[id]="."
      #data$diaphragm_E115[id]="."
      
    }
  }
  print("Lung Expression, the higher the more related")
  return(data)
}


gene_zscore<-function(data,lib){
  for(id in 1:dim(data)[1]){
    gid<-which(names(data)=="Gene.refGene"|names(data)=="GeneName")[1]
    index<-grep(data[id,gid],lib$gene)[1]
    if(length(index)>0){
      data$pLI[id]=lib$pLI[index]
      data$mis_z[id]=lib$mis_z[index]
      data$lof_z[id]=lib$lof_z[index]
    }else{
      data$mis_z[id]="."
      data$pLI[id]="."
      data$lof_z[id]="."
    }
  }
  return(data)
}




seperate_proband_info<-function(data){
  a<-data$proband
  if(grepl(pattern = "(",x = a,fixed = T)==F||grepl(pattern = ":",x = a,fixed = T)==F){ return(data)}
  probs<-mapply(a,FUN = function(x) unlist(strsplit(unlist(strsplit (x,split=":"))[1],split = "(",fixed = T))[1] )
  gts<-mapply(a,FUN = function(x) unlist(strsplit(unlist(strsplit (x,split=":"))[1],split = "(",fixed = T))[2] )
  ads<-mapply(a,FUN = function(x) unlist(strsplit (x,split=":"))[2])
  dps<-mapply(a,FUN = function(x) as.numeric(unlist(strsplit (x,split=":"))[3]))
  gqs<-mapply(a,FUN = function(x) unlist(strsplit (x,split=":"))[4])
  
  ab<-c();
  
  for(i in 1:length(gts)){
    gt=gts[i];
    
    alt=max(as.numeric(unlist(strsplit(gt,split = "/",fixed = T))),na.rm = T)
    
    ad=as.numeric(unlist(strsplit(ads[i],split = ","))[alt+1])
    ab<-c(ab,ad/dps[i])
    
  }
  max_ad<-mapply(ads,FUN = function(x) max(unlist(strsplit(x,split = ",")) ))
  data$proband<-probs
  data$GT<-gts
  data$AD<-ads
  data$Ind_DP<-dps
  data$AB<-ab
  data$Max_AD<-max_ad
  data$GQ<-gqs
  return (data)
}


addPhenotype2<-function(data,pheno){
  for(id in 1:length(data$proband)){
    key=data$proband[id];
    item<-which(pheno$ID==key)
    nms<-names(pheno)
    nms<-nms[which(nms!="")]
    if(length(item)>0){
      for(m in nms){
        data[id,as.character(m)]<-pheno[item,as.character(m)]
      }
    }
  }
  return (list(data=data,items=nms))
}

addPhenotype<-function(data,pheno){
  for(id in 1:length(data$proband)){
    key=unlist(strsplit(data$proband[id],split="_"))[2];
    item<-which(pheno$ID==key)
    nms<-names(pheno)
    nms<-nms[which(nms!="")]
    if(length(item)>0){
      for(m in nms){
        data[id,as.character(m)]<-pheno[item,as.character(m)]
      }
    }
  }
  return (list(data=data,items=nms))
}

output_candidates_v0<-function(denovo6,file2,items){
  if(missing(items)){items<-c()}
  denovo6$VarType="."
  denovo6$VarType[which(denovo6$MetaSVMprd=="D")]<-"DMIS"
  
  denovo6$VarType[which(denovo6$CADDphred>15 &denovo6$PP2.hdiv.prd=="D")]<-"DMIS"
  denovo6$VarType[which(denovo6$CADDphred>15 &denovo6$MetaSVMprd=="T")]<-"PDMIS"
  
  denovo6$VarType[which(denovo6$VarClass %in%c("stopgain","stoploss","frameshiftinsertion","frameshiftdeletion"))]<-"LGD"
  
  denovo6$VarType[which(denovo6$CADDphred<15 &denovo6$VarClass=="nonsynonymousSNV")]<-"MIS"
  
  aa<-mapply(denovo6$AAChange,FUN = function(x) unlist(strsplit(x,":"))[5])
  denovo6$AAC<-aa
  Nucleo<-mapply(denovo6$AAChange,FUN = function(x) unlist(strsplit(x,":"))[4])
  denovo6$Nucleotide<-Nucleo
  
  trans<-mapply(denovo6$AAChange,FUN = function(x) unlist(strsplit(x,":"))[2])
  denovo6$Transcript<-trans
  outsets<-c("proband","parents","CHROM","ID","POS","REF","ALT","FILTER","VarType","GeneName","VarClass","Transcript","Nucleotide","AAC","CADDphred","MetaSVMprd",
             "PP2.hdiv.prd","ExACfreq","1KGfreq","ESPfreq","75bp_mappability","AD","Ind_DP","GQ","IGV",
             "pLI","mis_z","LUNG_EXP","HEART_EXP","state","Gender","AGE","DIAG_AGE","pheno","group","population","AF","AC","AN",items)
  outset<-intersect(outsets,names(denovo6))
  write.csv(denovo6[,outset],
            file =file2,quote = T,row.names = F)
  
  return(denovo6)
}




output_candidates_v1<-function(denovo6,file2,items){
  if(missing(items)){items<-c()}
  denovo6$VarType="."
  denovo6$VarType[which(denovo6$MetaSVM_pred=="D")]<-"DMIS"
  
  denovo6$VarType[which(denovo6$CADD13_PHRED>15 &denovo6$Polyphen2_HVAR_pred=="D")]<-"DMIS"
  denovo6$VarType[which(denovo6$CADD13_PHRED>15 &denovo6$MetaSVM_pred=="T")]<-"PDMIS"
  
  denovo6$VarType[which(denovo6$ExonicFunc.refGene %in%c("stopgain","stoploss","frameshift_insertion","frameshift_deletion"))]<-"LGD"
  
  denovo6$VarType[which(denovo6$CADD13_PHRED<15 &denovo6$ExonicFunc.refGene=="nonsynonymous_SNV")]<-"MIS"
  
  aa<-mapply(denovo6$AAChange.refGene,FUN = function(x) unlist(strsplit(x,":"))[5])
  denovo6$AAC<-aa
  Nucleo<-mapply(denovo6$AAChange.refGene,FUN = function(x) unlist(strsplit(x,":"))[4])
  denovo6$Nucleotide<-Nucleo
  
  trans<-mapply(denovo6$AAChange.refGene,FUN = function(x) unlist(strsplit(x,":"))[2])
  denovo6$Transcript<-trans
  outsets<-c("proband","VarType","ExonicFunc.refGene","Gene.refGene","Transcript","Nucleotide","AAC","pLI","lof_z","mis_z","HEART_EXP","LUNG_EXP",
             "MetaSVM_pred", "Polyphen2_HDIV_pred","CADD13_PHRED",
             "ExAC_ALL","DIAG_AGE","pheno","SIFT_pred","MCAP","ExAC_AFR","ExAC_NFE","ExAC_AMR","ExAC_EAS","ExAC_OTH",
             "1000g2015aug_all","1000g2015aug_eas","1000g2015aug_eur",
             "AF","AC","AB","Ind_DP","QUAL","proband","state","Gender","AGE","CHROM","ID","POS","REF","ALT","FILTER","parents",items)
  outsets2<-c("proband","GeneName","VarClass","CADDphred","MetaSVMprd",
              "PP2.hdiv.prd","ExACfreq","1KGfreq","ESPfreq","Mappability","pheno","group","population")
  outsets<-c(outsets,outsets2)
  outset<-intersect(outsets,names(denovo6))
  write.csv(denovo6[,outset],
            file =file2,quote = T,row.names = F)
  
  return(denovo6)
}



filter_allfreq_local <- function(data,freq_avg,freq_max){
  
  data <- data[which(
    na.pass(as.numeric(data$ExAC_ALL)<= freq_avg)
    # &na.pass(as.numeric(data$ExAC_EUR)< freq_max)
  #  &na.pass(as.numeric(data$esp6500siv2_all)<= freq_avg)
    &as.numeric(data$ExAC_NFE)<= freq_max
    & as.numeric(data$gnomAD_exome_ALL)<= freq_max
    &as.numeric(data$gnomAD_exome_NFE)<= freq_max
    
    & as.numeric(data$`1000g2015aug_all`) <= freq_avg
    & as.numeric(data$`1000g2015aug_eur`) <= freq_avg
    
  ),]
  # data <- data[which( as.numeric(data$AC)< 25
  #                      &as.numeric(data$AB)>0.2
  #  ),]
  
  return(data)
}



COV_filter <- function(data){
  if(length(grep("gnmad_genome_",names(data)))>0){
    data<-data[which(as.numeric(data$gnmad_genome_10)>0.85),]
  }
  # if(length(grep("AC_PCGC",names(data)))>0 && length(grep("AC_xGen",names(data)))>0 && length(grep("AC_VCR",names(data)))>0){
  #   data <- data[which(2*data$AC_PCGC/data$AN_PCGC >0.8 & 2*data$AC_xGen/data$AN_xGen >0.8 & 2*(data$AC_VCR/data$AN_VCR) >0.8  &data$GnomAD_Genome_cov10>0.8 ),]
  # }
  # if(length(grep("AC_Control",names(data)))>0 && length(grep("AC_VU",names(data)))>0){
  #   data <- data[which(2*data$AC_Control/data$AN_Control >0.8 & 2*data$AC_VU/data$AN_VU >0.8 ),]
  # }
  return(data)
}



map_filter<-function(data){
  index<-grep("mappability",names(data),ignore.case = T)
  if(length(index)>0){ 
    data <- data[which(as.numeric(data[,index])==1|(data[,index]==".")),]
  }
  if(length(grep("genomicSuperDups",names(data),ignore.case = T))>0){
    index<-grep("Score",data$genomicSuperDups)
    
    dup_indexs<-index[which(as.numeric(unlist(lapply(index,FUN = function(x) unlist(strsplit(x = data$genomicSuperDups[x],split = ":|-|=|;|\\\\x3d|\\\\x3b"))[2])))>0.95)]
    if(length(dup_indexs)>0){
      data <- data[-dup_indexs,]
    }
  }
  return (data)
}



filter_case <- function(dat,total_case,f){
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  filter_rs<-c()
  print(dim(dat))
  max_AN<-max(dat$AN)
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("European",cnt))
  print(paste("European"))
  print(cnt)
  dat<-dat[which((dat$AN>0.9*max_AN & dat$`#CHROM`!="X")|dat$`#CHROM`=="X"),]
  dat<-dat[which((dat$AN>0.8*max_AN )),]
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("missing 0.1 ",cnt))
  print(paste("missingness <0.1"))
  print(cnt)
  
  #dat<-VQSR_filter(dat) 
  dat<-dat[which( dat$gnmad_genome_10!="." & (as.numeric(dat$gnmad_genome_15)>0.9 ) ),] # &  as.numeric(dat$gnomad_COV15)>0.75
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  print("90% cov10 and 75% cov15")
  filter_rs<-rbind(filter_rs,c("90% cov10 and 75% cov15 ",cnt))
  print(cnt)
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_case,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_case ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  print("Gene")
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  
  
  dat <- formatFreq_new(dat)
  print("Format freq")
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  print(cnt)
  
  dat <- filter_allfreq_local(dat,0.001,f) ## 
  # dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
  dat<-map_filter(dat)
#  dat<-dat[which( as.numeric(dat$BRAVO_AF)<f ),]
  # cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
  #        dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  # filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  # print("Bravo 10-4")
  # print(cnt)
  # 
  dat<-dat[which( as.numeric(dat$AF)< 100*f &  (as.numeric(dat$AC)/as.numeric(dat$AN)) < 100*f),]
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("cohorts 10-2 ",cnt))
  print("local cohort 10-2")
  print(cnt)
  
  #gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
 # dat$GQ<-gq
  rem<-which(as.numeric(dat$GQ_Ind)<99)
  if(length(rem)>0){
    PAH_CHD_dat<-dat[-rem,]  ## only for RGN
  }else{
    PAH_CHD_dat<-dat
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  print("gq 30")
  print(cnt)
  
  
  filter_rs<-rbind(filter_rs,c("gq 90 ",cnt))

  
 # PAH_CHD_dat<-ind_var(PAH_CHD_dat)
  # gseg<-unlist(lapply(PAH_CHD_dat$genomicSuperDups,FUN = function(x){if(x=="."){return (x);}else{ unlist(strsplit(unlist(strsplit(x,split = "="))[2],split = ";"))[1];} }))
  # rm<- which(gseg!="." & as.numeric(gseg) <0.95)
  # if(length(rm)>0){
  #   PAH_CHD_dat<-PAH_CHD_dat[-rm,]
  # }
  remove<- which( PAH_CHD_dat$AB< 0.3 |PAH_CHD_dat$AD< 8|PAH_CHD_dat$DP_ind<16)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  
  remove<- which( PAH_CHD_dat$Is_indel ==T  & PAH_CHD_dat$AB< 0.35)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  
  remove<- which( PAH_CHD_dat$Is_indel ==T  & PAH_CHD_dat$Is_homo==1 & PAH_CHD_dat$AB<1)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  remove<- which( PAH_CHD_dat$Is_indel ==T  & PAH_CHD_dat$Is_homo==0 & PAH_CHD_dat$AB > 0.9)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("AB 0.3 ",cnt))
  print("AB0.2")
  print(cnt)
  
  # PAH_CHD_dat$Is_indel=is_indel;
 # PAH_CHD_dat<-PAH_CHD_dat[which((PAH_CHD_dat$Is_indel==T)|  as.numeric(PAH_CHD_dat$GQ_Ind) >90 ),]
  # cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
  #        dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
  #        dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
  #        dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  # 
  # filter_rs<-rbind(filter_rs,c("GQ 90 ",cnt))
 
   # PAH_CHD_dat<-PAH_CHD_dat[which((PAH_CHD_dat$Is_indel==T)|  (PAH_CHD_dat$AB) >0.25 ),]
  
   # cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
   #        dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
   #        dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
   #        dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
   # 
   # filter_rs<-rbind(filter_rs,c("AB 0.3 SNV  ",cnt))
   # 
  # PAH_CHD_dat<-PAH_CHD_dat[which((PAH_CHD_dat$Is_indel==T)|  (PAH_CHD_dat$AD) >6 ),]
  PAH_CHD_dat<-PAH_CHD_dat[which(PAH_CHD_dat$FS <25),]
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("AD >6 SNV ",cnt))
#  PAH_CHD_dat$Is_indel<-is_indel
  indel_fail<-which(PAH_CHD_dat$FS>200|(PAH_CHD_dat$ReadPosRankSum!="." & as.numeric(PAH_CHD_dat$ReadPosRankSum) < -20) |PAH_CHD_dat$SOR > 10 )
  snv_fail<-c(which( (as.numeric(PAH_CHD_dat$MQ) < 35|
                                               as.numeric(PAH_CHD_dat$QD) < 2
                               
                                               | as.numeric(PAH_CHD_dat$SOR) > 6
                           
  )), which(  as.numeric(PAH_CHD_dat$ReadPosRankSum) < -2| as.numeric(PAH_CHD_dat$MQRankSum ) < -2 )) #| PAH_CHD_dat$SOR > 3
  # # snv_gq_fail<-c() #intersect(which(gq<70),which(is_indel==F))
  #  remove<- intersect(unique(snv_fail),which(PAH_CHD_dat$Is_indel==F)) #(c(indel_fail,snv_fail,snv_gq_fail))
  if(length(snv_fail)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-snv_fail,] ## only for RGN
  }
  
  
  # 
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("follow ",cnt))
  
  print(cnt )
  return(list(data=PAH_CHD_dat,filter=filter_rs))
}



filter_case_VQSR <- function(dat,total_case,f){
  dat<-dat[which(dat$ExonicFunc.refGene!="unknown"),]
  filter_rs<-c()
  print(dim(dat))
  max_AN<-max(dat$AN)
  
  cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  filter_rs<-rbind(filter_rs,c("European",cnt))
  print(paste("European"))
  print(cnt)
  dat<-dat[which((dat$AN>0.9*max_AN & dat$`#CHROM`!="X")|dat$`#CHROM`=="X"),]
  dat<-dat[which((dat$AN>0.8*max_AN )),]
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("missing 0.1 ",cnt))
  print(paste("missingness <0.1"))
  print(cnt)
  
  dat<-dat[grep("^exonic|^splic",dat$Func.refGene,ignore.case = T),]
  
  
  #dat<-VQSR_filter(dat) 
  dat<-dat[which( dat$gnmad_genome_10!="." & (as.numeric(dat$gnmad_genome_10)>0.9 ) ),] # &  as.numeric(dat$gnomad_COV15)>0.75
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
   filter_rs<-rbind(filter_rs,c("90% cov10  ",cnt))
  print(cnt)
  dat<-map_filter(dat)  ## mappability ==1
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("map ",cnt))
  print("map")
  print(c(dim(dat)/total_case,dim(dat[grep("^synony",dat$ExonicFunc.refGene),])/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])/total_case ))
  
  dat <-dat[grep("^MUC|^HLA",dat$Gene.refGene,invert = T),]
  print("Gene")
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("MUC|HLA ",cnt))
  
  
  dat <- formatFreq_new(dat)
  print("Format freq")
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  print(cnt)
  
  dat <- filter_allfreq_local(dat,0.01,f) ## 
  # dat <- dat[which(as.numeric(dat$gnomAD_genome_ALL)<0.0001 & as.numeric(dat$gnomAD_genome_NFE)<0.0001),] ## 
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("max pop 10-4 ",cnt))
  print("max pop 10-4")
  print(cnt)
 # dat<-map_filter(dat)
  #  dat<-dat[which( as.numeric(dat$BRAVO_AF)<f ),]
  # cnt<-c(dim( dat)[1]/total_case,dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1],
  #        dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  # filter_rs<-rbind(filter_rs,c("bravo 10-4 ",cnt))
  # print("Bravo 10-4")
  # print(cnt)
  # 
  dat<-dat[which( as.numeric(dat$AF)< 100*f &  (as.numeric(dat$AC)/as.numeric(dat$AN)) < 100*f),]
  cnt<-c(dim( dat)[1]/total_case,
         dim( dat[grep("^synony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim( dat[grep("^nonsynony", dat$ExonicFunc.refGene),])[1]/total_case,
         dim(dat[grep("non",dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("cohorts 10-2 ",cnt))
  print("local cohort 10-2")
  print(cnt)
  
  #gq<-unlist(lapply(dat$Genotype,FUN = function(x){unlist(strsplit(x,":"))[4]}))  ## For RGN
  # dat$GQ<-gq
  rem<-which(as.numeric(dat$GQ_Ind)< 90)
  if(length(rem)>0){
    PAH_CHD_dat<-dat[-rem,]  ## only for RGN
  }else{
    PAH_CHD_dat<-dat
  }
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  print("gq 70")
  print(cnt)
  PAH_CHD_dat<-PAH_CHD_dat[which(as.numeric(PAH_CHD_dat$MQ) > 40),]
  
  filter_rs<-rbind(filter_rs,c("gq 70 ",cnt))
  
  
  # PAH_CHD_dat<-ind_var(PAH_CHD_dat)
  # gseg<-unlist(lapply(PAH_CHD_dat$genomicSuperDups,FUN = function(x){if(x=="."){return (x);}else{ unlist(strsplit(unlist(strsplit(x,split = "="))[2],split = ";"))[1];} }))
  # rm<- which(gseg!="." & as.numeric(gseg) <0.95)
  # if(length(rm)>0){
  #   PAH_CHD_dat<-PAH_CHD_dat[-rm,]
  # }
  remove<- which( PAH_CHD_dat$AD< 5)  #PAH_CHD_dat$DP_ind<10 | PAH_CHD_dat$AB< 0.3 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  
  remove<- which( PAH_CHD_dat$Is_indel ==T  & PAH_CHD_dat$Is_homo==1 & PAH_CHD_dat$AB<1)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  remove<- which( PAH_CHD_dat$Is_indel ==T  & PAH_CHD_dat$Is_homo==0 & PAH_CHD_dat$AB > 0.9)  #PAH_CHD_dat$DP_ind<10 |
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  PAH_CHD_dat<- PAH_CHD_dat[which( as.numeric(PAH_CHD_dat$AB) > 0.25),]
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("AD 5 ",cnt))
  print("AB0.25")
  print(cnt)
  
  PAH_CHD_dat<- PAH_CHD_dat[which( as.numeric(PAH_CHD_dat$FS) < 20),]
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("FS 20 ",cnt))
  print("AB0.2")
  print(cnt)
  # rm<-intersect(which( PAH_CHD_dat$VQSLOD < -2 ),grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T))
  # if(length(rm)>0){
  #   PAH_CHD_dat<-PAH_CHD_dat[-rm,]
  # }
  # -rbind(filter_rs,c("follow ",cnt))
  # 
  
  cnt<-c(dim( PAH_CHD_dat)[1]/total_case,
         dim( PAH_CHD_dat[grep("^synony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim( PAH_CHD_dat[grep("^nonsynony", PAH_CHD_dat$ExonicFunc.refGene),])[1]/total_case,
         dim(PAH_CHD_dat[grep("non",PAH_CHD_dat$ExonicFunc.refGene,invert = T),])[1]/total_case) 
  
  filter_rs<-rbind(filter_rs,c("VQSLOD -5 ",cnt))
  print(cnt )
  
  for(lod in seq(-10,14,0.5)){
    tm<- PAH_CHD_dat[which(as.numeric(PAH_CHD_dat$VQSLOD) > lod),];
    cnt<-c(
      qcase(tm[grep("^synony",tm$ExonicFunc.refGene),c("proband","Gene.refGene")])/total_case,
      qcase(tm[grep("^nonsynony",tm$ExonicFunc.refGene),c("proband","Gene.refGene")])/total_case,
      qcase(tm[grep("non",tm$ExonicFunc.refGene,invert = T),c("proband","Gene.refGene")])/total_case)
    filter_rs<-rbind(filter_rs,c(paste("lod ",lod,": "),cnt))
  }
  
  PAH_CHD_dat<- PAH_CHD_dat[which( as.numeric(PAH_CHD_dat$VQSLOD) > -4),]
  rm<- which( PAH_CHD_dat$VQSLOD < 2 & PAH_CHD_dat$Is_indel == F)
  if(length(rm) >0){
    PAH_CHD_dat<- PAH_CHD_dat[-rm,]
  }
  
  snv_fail<-c(which( (as.numeric(PAH_CHD_dat$QD) < 1| as.numeric(PAH_CHD_dat$SOR) > 3)), 
              which(as.numeric(PAH_CHD_dat$ReadPosRankSum) < -2| as.numeric(PAH_CHD_dat$MQRankSum ) < -2 )) #| PAH_CHD_dat$SOR > 3
  # # snv_gq_fail<-c() #intersect(which(gq<70),which(is_indel==F))
  remove<- intersect(unique(snv_fail),which(PAH_CHD_dat$Is_indel==F)) #(c(indel_fail,snv_fail,snv_gq_fail))
  if(length(remove)>0){
    PAH_CHD_dat<-PAH_CHD_dat[-remove,] ## only for RGN
  }
  
#  PAH_CHD_dat<- PAH_CHD_dat[which( PAH_CHD_dat$VQSLOD > 4),]
  return(list(data=PAH_CHD_dat,filter=filter_rs))
}

filter_gnomad<-function(dat,f){
  total_control=7509
  gfilter_rs<-c()
  #chd_dat <- read.csv("PAH/PAH-CHD/CHD/CHD.WES.253known.csv",header = 1,stringsAsFactors = F,check.names = F,comment.char = "")
  subdat<-dat[which(dat$ExonicFunc.refGene!="unknown" & as.numeric(dat$AC_NFE)>0),] #,"RF" dat$FILTER%in%c("PASS") & 
  subdat<-subdat[grep("InbreedingCoeff",subdat$FILTER,invert = T),]

  cnt<- c(
           sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
           sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
           sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  
  gfilter_rs<-rbind(gfilter_rs,c("PASS",cnt))
  #  if(subdat$`#CHROM`!="X"){
  max_gan<-max(subdat$AN)
  index<-which((subdat$`#CHROM`!="X" & subdat$AN/max_gan >0.75)|subdat$`#CHROM`=="X")
  subdat<-subdat[index,]
  # subdat<-subdat[which(subdat$AN_NFE/max_gan >0.75),]
  
  # }
  cnt<- c(
    sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  
  gfilter_rs<-rbind(gfilter_rs,c("PASS",cnt))
  subdat<-subdat[which(as.numeric(subdat$DP_10) >=0.9),]
  subdat<-map_filter(subdat)
  cnt<- c(
    sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  
  gfilter_rs<-rbind(gfilter_rs,c("MAP",cnt))
  
  subdat <-subdat[grep("^MUC|^HLA",subdat$Gene.refGene,invert = T),]
  cnt<-c(
    sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  
  
  gfilter_rs<-rbind(gfilter_rs,c("MUC|HLA",cnt))
  
#  subdat<-map_filter(subdat)
  cnt<-c(
     sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
     sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
     sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
   
  gfilter_rs<-rbind(gfilter_rs,c("MAP",cnt))
  subdat <- formatFreq_new(subdat)
  
  subdat <- filter_allfreq_local(subdat,0.01,f)  ## Max gnomad exome <10^-4 ## Max gnomad exome <10^-4
  cnt<-c(
    sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  
  
  gfilter_rs<-rbind(gfilter_rs,c("MAX pop AF",cnt))
  
  
  
  subdat <- subdat[which(as.numeric(subdat$AF)< 100*f & as.numeric(subdat$AF_NFE)< 100*f  ),] 
  cnt<-c(
    sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
    sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  
  gfilter_rs<-rbind(gfilter_rs,c("local AF",cnt))
  print(cnt)
 
 
  # subdat <- subdat[which(as.numeric(subdat$FS)< 40  ),] 
  # cnt<-c(
  #   sum(subdat$AC_NFE[grep("^synony",subdat$ExonicFunc.refGene)])/total_control,
  #   sum(subdat$AC_NFE[grep("^nonsynony",subdat$ExonicFunc.refGene)])/total_control,
  #   sum(subdat$AC_NFE[grep("non",subdat$ExonicFunc.refGene,invert = T)])/total_control)
  # 
  # gfilter_rs<-rbind(gfilter_rs,c("FS ",cnt))
  # print(cnt)
  
  for(lod in seq(-10,4,0.5)){
    tm<- subdat[which(as.numeric(subdat$VQSLOD) > lod),];
    cnt<-c(
      sum(tm$AC_NFE[grep("^synony",tm$ExonicFunc.refGene)])/total_control,
      sum(tm$AC_NFE[grep("^nonsynony",tm$ExonicFunc.refGene)])/total_control,
      sum(tm$AC_NFE[grep("non",tm$ExonicFunc.refGene,invert = T)])/total_control)
    gfilter_rs<-rbind(gfilter_rs,c(paste("lod ",lod,": "),cnt))
  }
  subdat<-subdat[which(as.numeric(subdat$VQSLOD) > -11),];
  rm<-intersect(which(as.numeric(subdat$VQSLOD) < -5),grep("non",subdat$ExonicFunc.refGene,invert = T) );
  if(length(rm) >0){
     subdat<-subdat[-rm,]
   }
  return (list(data=subdat,filter=gfilter_rs))
}

