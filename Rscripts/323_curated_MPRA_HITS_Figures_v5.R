#

suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
library("ggridges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL


ADD_INTERVAL_DE_Main_VARS = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
 

  #### CURATED TABLE ----
  
  CURATED_TABLE<-as.data.frame(fread(file=opt$CURATED_TABLE, sep =",", header =T), stringsAsFactors=F)
  
  colnames(CURATED_TABLE)<-gsub(" +","_",colnames(CURATED_TABLE))
  # colnames(CURATED_TABLE)[4]<-"Candidate_effector"
  
  #CURATED_TABLE<-gsub(" +","_",CURATED_TABLE)
  
  
  # CURATED_TABLE$Mechanistic_Class[which(CURATED_TABLE$Mechanistic_Class == "")]<-"No regulation"
  # CURATED_TABLE$Candidate_effector[which(CURATED_TABLE$Candidate_effector == "")]<-"ABSENT"
  # CURATED_TABLE$OpenTargets_QTL[which(CURATED_TABLE$OpenTargets_QTL == "")]<-"NOVEL"
  
  
  
  cat("CURATED_TABLE_0:\n")
  cat(str(CURATED_TABLE))
  cat("\n")
  cat(str(unique(CURATED_TABLE$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(CURATED_TABLE$Candidate_effector))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CURATED_TABLE$Candidate_effector)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(CURATED_TABLE$Mechanistic_Class))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CURATED_TABLE$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(CURATED_TABLE$OpenTargets_QTL))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(CURATED_TABLE$OpenTargets_QTL)))))
  cat("\n")
  
  ### Read Results_INTERVAL_LM----
  
  
  
  Results_INTERVAL_LM<-as.data.frame(fread(file=opt$Results_INTERVAL_LM, sep="\t", header=T) , stringsAsFactors=T)
  
  
  cat("Results_INTERVAL_LM\n")
  cat(str(Results_INTERVAL_LM))
  cat("\n")
  cat(str(unique(Results_INTERVAL_LM$VAR)))
  cat("\n")
  
  
  ####MASTER LOOP----
  
  Condition_DEBUG <- 0
  
  VARS<-unique(Results_INTERVAL_LM$VAR)
  
  list_DEF<-list()
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("-------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    Results_INTERVAL_LM_sel<-Results_INTERVAL_LM[which(Results_INTERVAL_LM$VAR == VAR_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Results_INTERVAL_LM_sel_\n")
      cat(str(Results_INTERVAL_LM_sel))
      cat("\n")
      cat(str(unique(Results_INTERVAL_LM_sel$ensembl_gene_id)))
      cat("\n")
    }
    
    ENSG_array<-unique(c(Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$CIS_gene_minuslogpvalue)],
                         Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$Block_PCHiC_minuslogpvalue)]))
    
    
    
    if(length(ENSG_array) >0)
    {
      # cat("ENSG_array_\n")
      # cat(str(ENSG_array))
      # cat("\n")
      
      
      Results_INTERVAL_LM_sel_ENSG_sel<-Results_INTERVAL_LM_sel[which(Results_INTERVAL_LM_sel$ensembl_gene_id%in%ENSG_array),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_INTERVAL_LM_sel_ENSG_sel_\n")
        cat(str(Results_INTERVAL_LM_sel_ENSG_sel))
        cat("\n")
        cat(str(unique(Results_INTERVAL_LM_sel_ENSG_sel$ensembl_gene_id)))
        cat("\n")
      }
      
      CIS_gene_INTERVAL<-Results_INTERVAL_LM_sel_ENSG_sel[!is.na(Results_INTERVAL_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                        !is.na(Results_INTERVAL_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
      
      if(Condition_DEBUG == 1)
      {
        cat("CIS_gene_\n")
        cat(str(CIS_gene_INTERVAL))
        cat("\n")
        cat(str(unique(CIS_gene_INTERVAL$ensembl_gene_id)))
        cat("\n")
      }
      
      Block_genes_INTERVAL<-Results_INTERVAL_LM_sel_ENSG_sel[is.na(Results_INTERVAL_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                           !is.na(Results_INTERVAL_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Block_genes_\n")
        cat(str(Block_genes_INTERVAL))
        cat("\n")
        cat(str(unique(Block_genes_INTERVAL$ensembl_gene_id)))
        cat("\n")
      }
      
      #"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS","coefficient_Genotypes_specific_CELL_COUNTS",
      
      DE_gene_df<-band.genes<-data.frame(matrix(ncol=4,nrow=length(ENSG_array), 
                                                dimnames=list(NULL, c("VAR", "ensembl_gene_id",
                                                                      "Significance",
                                                                      "RNASeq_source"))),
                                         stringsAsFactors = F)
      
      DE_gene_df$VAR<-VAR_sel
      DE_gene_df$ensembl_gene_id<-ENSG_array
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_gene_df_0\n")
        cat(str(DE_gene_df))
        cat("\n")
        cat(str(unique(DE_gene_df$ensembl_gene_id)))
        cat("\n")
      }
      
      if(dim(CIS_gene_INTERVAL)[1] >0)
      {
        Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                             which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                             which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                             which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"))]
        colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
        colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
        
        if(Condition_DEBUG == 1)
        {
          cat("Block_genes_INTERVAL_subset_0\n")
          cat(str(Block_genes_INTERVAL_subset))
          cat("\n")
          
        }
        
        CIS_gene_INTERVAL_subset<-CIS_gene_INTERVAL[,c(which(colnames(CIS_gene_INTERVAL) == "ensembl_gene_id"),
                                                       which(colnames(CIS_gene_INTERVAL) == "HGNC"),
                                                       which(colnames(CIS_gene_INTERVAL) == "CIS_gene_minuslogpvalue"),
                                                       which(colnames(CIS_gene_INTERVAL) == "CIS_Beta"))]
        colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
        colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
        
        if(Condition_DEBUG == 1)
        {
          cat("CIS_gene_INTERVAL_subset_0\n")
          cat(str(CIS_gene_INTERVAL_subset))
          cat("\n")
          
        }
        
        
        Bind_CIS_Block_INTERVAL<-rbind(CIS_gene_INTERVAL_subset,
                                       Block_genes_INTERVAL_subset)
        
        if(Condition_DEBUG == 1)
        {
          cat("Bind_CIS_Block_INTERVAL_0\n")
          cat(str(Bind_CIS_Block_INTERVAL))
          cat("\n")
          
        }
        
        DE_gene_df<-merge(DE_gene_df,
                          Bind_CIS_Block_INTERVAL,
                          by=c("ensembl_gene_id"),
                          all=T)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_gene_df_1\n")
          cat(str(DE_gene_df))
          cat("\n")
          cat(str(unique(DE_gene_df$ensembl_gene_id)))
          cat("\n")
        }
        
        # quit(status = 1)
        
      }else{
        
        Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                             which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                             which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                             which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"))]
        colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
        colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
        
        if(Condition_DEBUG == 1)
        {
          cat("Block_genes_INTERVAL_subset_0\n")
          cat(str(Block_genes_INTERVAL_subset))
          cat("\n")
          
        }
        
        DE_gene_df<-merge(DE_gene_df,
                          Block_genes_INTERVAL_subset,
                          by=c("ensembl_gene_id"),
                          all=T)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_gene_df_1\n")
          cat(str(DE_gene_df))
          cat("\n")
          cat(str(unique(DE_gene_df$ensembl_gene_id)))
          cat("\n")
        }
        
      }# dim(CIS_gene_INTERVAL)[1] >0
      
      
      #### Significance ----
      
      DE_gene_df$Significance[which(DE_gene_df$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3)]<-"YES"
      DE_gene_df$Significance[which(DE_gene_df$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS < 1.3)]<-"NO"
      
      DE_gene_df$Significance<-factor(DE_gene_df$Significance,
                                      levels=c("NO","YES"),
                                      ordered=T)
      if(Condition_DEBUG == 1)
      {
        cat("DE_gene_df$Significance\n")
        cat(sprintf(as.character(names(summary(DE_gene_df$Significance)))))
        cat("\n")
        cat(sprintf(as.character(summary(DE_gene_df$Significance))))
        cat("\n")
      }
      
      #### RNASeq_source ----
      
      DE_gene_df$RNASeq_source<-"Whole blood"
      
      
      DE_gene_df$RNASeq_source<-factor(DE_gene_df$RNASeq_source,
                                       levels=c("Whole blood"),
                                       ordered=T)
      if(Condition_DEBUG == 1)
      {
        cat("DE_gene_df$RNASeq_source\n")
        cat(sprintf(as.character(names(summary(DE_gene_df$RNASeq_source)))))
        cat("\n")
        cat(sprintf(as.character(summary(DE_gene_df$RNASeq_source))))
        cat("\n")
      }
      
      DE_gene_df_SIG<-DE_gene_df[which(DE_gene_df$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
      
      if(dim(DE_gene_df_SIG)[1] >0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("DE_gene_df_SIG_1\n")
          cat(str(DE_gene_df_SIG))
          cat("\n")
          cat(str(unique(DE_gene_df_SIG$ensembl_gene_id)))
          cat("\n")
        }
        
        DE_gene_df_SIG.dt<-data.table(DE_gene_df_SIG, key="VAR")
        
        DE_gene_df_SIG_collapsed<-as.data.frame(DE_gene_df_SIG.dt[,.(Whole_blood_DE_HGNC_string=paste(HGNC, collapse=";")),
                                                                  by=key(DE_gene_df_SIG.dt)], stringsAsFactors=F)
        
        if(Condition_DEBUG == 1)
        {
          cat("DE_gene_df_SIG_collapsed_1\n")
          cat(str(DE_gene_df_SIG_collapsed))
          cat("\n")
          
        }
        
        list_DEF[[i]]<-DE_gene_df_SIG_collapsed
        # ############################
        # quit(status = 1)
        
      }#dim(DE_gene_df_SIG)[1] >0
    }#length(ENSG_array) >0
  }#i in 1:length(VARS)
  
  Condition_DEBUG <- 1
  
  if(length(list_DEF) >0)
  {
    
    DE_gene_df_SIG_collapsed_DEF = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
    
    if(Condition_DEBUG == 1)
    {
      cat("DE_gene_df_SIG_collapsed_DEF_0\n")
      cat(str(DE_gene_df_SIG_collapsed_DEF))
      cat("\n")
    }
    
    CURATED_TABLE<-merge(CURATED_TABLE,
                         DE_gene_df_SIG_collapsed_DEF,
                         by="VAR",
                         all.x=T)
    
    if(Condition_DEBUG == 1)
    {
      cat("CURATED_TABLE_1\n")
      cat(str(CURATED_TABLE))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(CURATED_TABLE$Whole_blood_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(CURATED_TABLE$Whole_blood_DE_HGNC_string)))))
      cat("\n")
    }
    
    
    path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
   
    if (file.exists(path7)){
      
      
      
      
    } else {
      dir.create(file.path(path7))
      
    }
    
    if(Condition_DEBUG == 1)
    {
      cat("path7:\t")
      cat(sprintf(as.character(path7)))
      cat("\n")
    }
    
    setwd(path7)
    
    unlink("ZZZ_in_process.tsv")
    
    write.table(CURATED_TABLE, sep="\t", quote=F, row.names = F,file="ZZZ_in_process.tsv")
    
  }# length(list_DEF) >0
  
}

ADD_INTERVAL_DTU_Main_VARS = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  
  #### CURATED TABLE ----
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
 
  setwd(path7)
  
  ZZZ_table<-as.data.frame(fread(file="ZZZ_in_process.tsv", sep ="\t", header =T), stringsAsFactors=F)
  
 
  
  cat("ZZZ_table_0:\n")
  cat(str(ZZZ_table))
  cat("\n")
  cat(str(unique(ZZZ_table$VAR)))
  cat("\n")
 
  
  ### Read Results_LogLM----
  
  
  
  Results_LogLM<-as.data.frame(fread(file=opt$Results_LogLM, sep="\t", header=T) , stringsAsFactors=T)
  
  
  cat("Results_LogLM\n")
  cat(str(Results_LogLM))
  cat("\n")
  cat(str(unique(Results_LogLM$VAR)))
  cat("\n")
  
  
  ####MASTER LOOP----
  
  Condition_DEBUG <- 0
  
  VARS<-unique(Results_LogLM$VAR)
  
  list_DEF<-list()
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("-------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    Results_LogLM_sel<-Results_LogLM[which(Results_LogLM$VAR == VAR_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Results_LogLM_sel_\n")
      cat(str(Results_LogLM_sel))
      cat("\n")
      cat(str(unique(Results_LogLM_sel$ensembl_gene_id)))
      cat("\n")
    }
    
    ENSG_array<-unique(c(Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$CIS_gene_minuslogpvalue)],
                         Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$Block_PCHiC_minuslogpvalue)]))
    
    
    
    ENSG_array<-unique(c(Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$CIS_gene_minuslogpvalue)],
                         Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$Block_PCHiC_minuslogpvalue)]))
    
    
    
    if(length(ENSG_array) >0)
    {
      cat("ENSG_array_\n")
      cat(str(ENSG_array))
      cat("\n")
      
      
      list_per_gene<-list()
      for(z in 1:length(ENSG_array))
      {
        ENSG_array_sel<-ENSG_array[z]
        
        
        Results_LogLM_sel_ENSG_sel<-Results_LogLM_sel[which(Results_LogLM_sel$ensembl_gene_id%in%ENSG_array_sel),]
        
        HGNC_sel<-unique(Results_LogLM_sel_ENSG_sel$HGNC)
        
        cat("------->\t")
        cat(sprintf(as.character(HGNC_sel)))
        cat("\t")
        cat(sprintf(as.character(ENSG_array_sel)))
        cat("\n")
        
        
        
        if(Condition_DEBUG == 1)
        {
          cat("Results_LogLM_sel_ENSG_sel_\n")
          cat(str(Results_LogLM_sel_ENSG_sel))
          cat("\n")
          cat(str(unique(Results_LogLM_sel_ENSG_sel$ensembl_gene_id)))
          cat("\n")
        }
        
        CIS_gene_INTERVAL<-Results_LogLM_sel_ENSG_sel[!is.na(Results_LogLM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                        !is.na(Results_LogLM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
        
        if(Condition_DEBUG == 1)
        {
          cat("CIS_gene_\n")
          cat(str(CIS_gene_INTERVAL))
          cat("\n")
          cat(str(unique(CIS_gene_INTERVAL$ensembl_gene_id)))
          cat("\n")
        }
        
        Block_genes_INTERVAL<-Results_LogLM_sel_ENSG_sel[is.na(Results_LogLM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                           !is.na(Results_LogLM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
        
        if(Condition_DEBUG == 1)
        {
          cat("Block_genes_\n")
          cat(str(Block_genes_INTERVAL))
          cat("\n")
          cat(str(unique(Block_genes_INTERVAL$ensembl_gene_id)))
          cat("\n")
        }
        
        #"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS","coefficient_Genotypes_specific_CELL_COUNTS",
        
        df_aggregate<-data.frame()
        
        if(dim(CIS_gene_INTERVAL)[1] >0)
        {
          Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "transcript_id"),
                                                               which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                               which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                               which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                               which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"))]
          colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
          colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
          
          if(Condition_DEBUG == 1)
          {
            cat("Block_genes_INTERVAL_subset_0\n")
            cat(str(Block_genes_INTERVAL_subset))
            cat("\n")
            
          }
          
          CIS_gene_INTERVAL_subset<-CIS_gene_INTERVAL[,c(which(colnames(CIS_gene_INTERVAL) == "transcript_id"),
                                                         which(colnames(CIS_gene_INTERVAL) == "ensembl_gene_id"),
                                                         which(colnames(CIS_gene_INTERVAL) == "HGNC"),
                                                         which(colnames(CIS_gene_INTERVAL) == "CIS_gene_minuslogpvalue"),
                                                         which(colnames(CIS_gene_INTERVAL) == "CIS_Beta"))]
          colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
          colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
          
          if(Condition_DEBUG == 1)
          {
            cat("CIS_gene_INTERVAL_subset_0\n")
            cat(str(CIS_gene_INTERVAL_subset))
            cat("\n")
            
          }
          
          
          Bind_CIS_Block_INTERVAL<-rbind(CIS_gene_INTERVAL_subset,
                                         Block_genes_INTERVAL_subset)
          
          if(Condition_DEBUG == 1)
          {
            cat("Bind_CIS_Block_INTERVAL_0\n")
            cat(str(Bind_CIS_Block_INTERVAL))
            cat("\n")
            
          }
          
          df_aggregate<-rbind(df_aggregate,Bind_CIS_Block_INTERVAL)
          
        }else{
          
          Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "transcript_id"),
                                                               which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                               which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                               which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                               which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"))]
          colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
          colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
          
          if(Condition_DEBUG == 1)
          {
            cat("Block_genes_INTERVAL_subset_0\n")
            cat(str(Block_genes_INTERVAL_subset))
            cat("\n")
            
          }
          
          df_aggregate<-rbind(df_aggregate,Block_genes_INTERVAL_subset)
          
        }# dim(CIS_gene_INTERVAL)[1] >0
        
        if(Condition_DEBUG == 1)
        {
          cat("df_aggregate_0\n")
          cat(str(df_aggregate))
          cat("\n")
          
        }
        
        rescue_transcript_INTERVAL<-Results_LogLM_sel_ENSG_sel[-which(Results_LogLM_sel_ENSG_sel$transcript_id%in%df_aggregate$transcript_id),c(which(colnames(Results_LogLM_sel_ENSG_sel) == "transcript_id"),
                                                                                                                                                which(colnames(Results_LogLM_sel_ENSG_sel) == "ensembl_gene_id"),
                                                                                                                                                which(colnames(Results_LogLM_sel_ENSG_sel) == "HGNC"),
                                                                                                                                                which(colnames(Results_LogLM_sel_ENSG_sel) == "Block_PCHiC_minuslogpvalue"),
                                                                                                                                                which(colnames(Results_LogLM_sel_ENSG_sel) == "Block_PCHiC_Beta"))]
        colnames(rescue_transcript_INTERVAL)[which(colnames(rescue_transcript_INTERVAL)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
        colnames(rescue_transcript_INTERVAL)[which(colnames(rescue_transcript_INTERVAL)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
        
        if(Condition_DEBUG == 1)
        {
          cat("rescue_transcript_INTERVAL_0\n")
          cat(str(rescue_transcript_INTERVAL))
          cat("\n")
          
        }
        
        
        
        df_aggregate<-rbind(df_aggregate,rescue_transcript_INTERVAL)
        
        
        if(Condition_DEBUG == 1)
        {
          cat("df_aggregate_1\n")
          cat(str(df_aggregate))
          cat("\n")
          
        }
        
        df_aggregate$VAR<-VAR_sel
        df_aggregate$Significance<-NA
        df_aggregate$RNASeq_source<-NA
        
        if(Condition_DEBUG == 1)
        {
          cat("df_aggregate_2\n")
          cat(str(df_aggregate))
          cat("\n")
          
        }
        
        #### Significance ----
        
        df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3)]<-"YES"
        df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS < 1.3)]<-"NO"
        
        df_aggregate$Significance<-factor(df_aggregate$Significance,
                                          levels=c("NO","YES"),
                                          ordered=T)
        if(Condition_DEBUG == 1)
        {
          cat("df_aggregate$Significance\n")
          cat(sprintf(as.character(names(summary(df_aggregate$Significance)))))
          cat("\n")
          cat(sprintf(as.character(summary(df_aggregate$Significance))))
          cat("\n")
        }
        
        #### RNASeq_source ----
        
        df_aggregate$RNASeq_source<-"Whole blood"
        
        
        df_aggregate$RNASeq_source<-factor(df_aggregate$RNASeq_source,
                                           levels=c("Whole blood"),
                                           ordered=T)
        if(Condition_DEBUG == 1)
        {
          cat("df_aggregate_3\n")
          cat(str(df_aggregate))
          cat("\n")
          
        }
        
        df_aggregate.dt<-data.table(df_aggregate, key="VAR")
        
        
        df_aggregate_MAX<-as.data.frame(df_aggregate.dt[,.SD[which.max(ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS)],by=key(df_aggregate.dt)], stringsAsFactors=F)
        
        if(Condition_DEBUG == 1)
        {
          cat("df_aggregate_MAX_0\n")
          cat(str(df_aggregate_MAX))
          cat("\n")
          
        }
      
        df_aggregate_MAX_SIG<-df_aggregate_MAX[which(df_aggregate_MAX$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
        
        if(dim(df_aggregate_MAX_SIG)[1] >0)
        {
          if(Condition_DEBUG == 1)
          {
            cat("df_aggregate_MAX_SIG_0\n")
            cat(str(df_aggregate_MAX_SIG))
            cat("\n")
            
          }
          
          list_per_gene[[z]]<-df_aggregate_MAX_SIG
        }#dim(df_aggregate_MAX_SIG)[1] >0
      }#z in 1:length(ENSG_array)
      
      if(length(list_per_gene) >0)
      {
        DTU_genes_per_VAR = unique(as.data.frame(data.table::rbindlist(list_per_gene, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("DTU_genes_per_VAR_0\n")
          cat(str(DTU_genes_per_VAR))
          cat("\n")
        }
        
        DTU_genes_per_VAR.dt<-data.table(DTU_genes_per_VAR, key="VAR")
        
        DTU_genes_per_VAR_collapsed<-DTU_genes_per_VAR.dt[,.(Whole_blood_DTU_HGNC_string=paste(HGNC, collapse=";")),by=key(DTU_genes_per_VAR.dt)]
        
        if(Condition_DEBUG == 1)
        {
          cat("DTU_genes_per_VAR_collapsed_0\n")
          cat(str(DTU_genes_per_VAR_collapsed))
          cat("\n")
        }
        
        list_DEF[[i]]<-DTU_genes_per_VAR_collapsed
        
        # quit(status = 1)
      }#length(list_per_gene) >0
    }#length(ENSG_array) >0
  }#i in 1:length(VARS)
  
  Condition_DEBUG <- 1
  
  if(length(list_DEF) >0)
  {
    
    DTU_gene_df_SIG_collapsed_DEF = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
    
    if(Condition_DEBUG == 1)
    {
      cat("DTU_gene_df_SIG_collapsed_DEF_0\n")
      cat(str(DTU_gene_df_SIG_collapsed_DEF))
      cat("\n")
    }
    
    ZZZ_table<-merge(ZZZ_table,
                         DTU_gene_df_SIG_collapsed_DEF,
                         by="VAR",
                         all.x=T)
    
    if(Condition_DEBUG == 1)
    {
      cat("ZZZ_table_1\n")
      cat(str(ZZZ_table))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Whole_blood_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Whole_blood_DE_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Whole_blood_DTU_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Whole_blood_DTU_HGNC_string)))))
      cat("\n")
      
      check<-ZZZ_table[(is.na(ZZZ_table$Whole_blood_DE_HGNC_string) &
                         is.na(ZZZ_table$Whole_blood_DTU_HGNC_string)),]
      
      cat("check_1\n")
      cat(str(check))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(check$Whole_blood_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(check$Whole_blood_DE_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(check$Whole_blood_DTU_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(check$Whole_blood_DTU_HGNC_string)))))
      cat("\n")
      
      
    }
    
    
    path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
    
    if (file.exists(path7)){
      
      
      
      
    } else {
      dir.create(file.path(path7))
      
    }
    
    if(Condition_DEBUG == 1)
    {
      cat("path7:\t")
      cat(sprintf(as.character(path7)))
      cat("\n")
    }
    
    setwd(path7)
    
    write.table(ZZZ_table, sep="\t", quote=F, row.names = F,file="ZZZ_in_process.tsv")
    
  }# length(list_DEF) >0
  
}

ADD_BP_DE_Main_VARS = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  
  #### CURATED TABLE ----
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  setwd(path7)
  
  ZZZ_table<-as.data.frame(fread(file="ZZZ_in_process.tsv", sep ="\t", header =T), stringsAsFactors=F)
  
  
  
  cat("ZZZ_table_0:\n")
  cat(str(ZZZ_table))
  cat("\n")
  cat(str(unique(ZZZ_table$VAR)))
  cat("\n")
  
  
  ### Read Results_BP_LM----
  
  
  
  Results_BP_LM<-as.data.frame(fread(file=opt$Results_BP_LM, sep="\t", header=T) , stringsAsFactors=T)
  Results_BP_LM$Cell_Type<-factor(Results_BP_LM$Cell_Type, levels=c("Monocyte","Neutrophil","Tcell"), ordered=T)
  
  
  cat("Results_BP_LM\n")
  cat(str(Results_BP_LM))
  cat("\n")
  cat(str(unique(Results_BP_LM$VAR)))
  cat("\n")
  
  
  ####MASTER LOOP----
  
  Condition_DEBUG <- 0
  
  VARS<-unique(Results_BP_LM$VAR)
  
  list_DEF<-list()
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    cat("-------------------->\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    Results_BP_LM_sel<-Results_BP_LM[which(Results_BP_LM$VAR == VAR_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Results_BP_LM_sel_\n")
      cat(str(Results_BP_LM_sel))
      cat("\n")
      cat(str(unique(Results_BP_LM_sel$ensembl_gene_id)))
      cat("\n")
    }
    
    if(dim(Results_BP_LM_sel)[1] >0)
    {
      
      
      ENSG_array<-unique(c(Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$CIS_gene_minuslogpvalue)],
                           Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$Block_PCHiC_minuslogpvalue)]))
      
      
      
      ENSG_array<-unique(c(Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$CIS_gene_minuslogpvalue)],
                           Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$Block_PCHiC_minuslogpvalue)]))
      
      
      
      if(length(ENSG_array) >0)
      {
        cat("ENSG_array_\n")
        cat(str(ENSG_array))
        cat("\n")
        
        
        list_per_gene<-list()
        for(z in 1:length(ENSG_array))
        {
          ENSG_array_sel<-ENSG_array[z]
          
          
          Results_BP_LM_sel_ENSG_sel<-Results_BP_LM_sel[which(Results_BP_LM_sel$ensembl_gene_id%in%ENSG_array_sel),]
          
          HGNC_sel<-unique(Results_BP_LM_sel_ENSG_sel$HGNC)
          
          cat("------->\t")
          cat(sprintf(as.character(HGNC_sel)))
          cat("\t")
          cat(sprintf(as.character(ENSG_array_sel)))
          cat("\n")
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_BP_LM_sel_ENSG_sel_\n")
            cat(str(Results_BP_LM_sel_ENSG_sel))
            cat("\n")
            cat(str(unique(Results_BP_LM_sel_ENSG_sel$ensembl_gene_id)))
            cat("\n")
          }
          
          CIS_gene_INTERVAL<-Results_BP_LM_sel_ENSG_sel[!is.na(Results_BP_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                          !is.na(Results_BP_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
          
          if(Condition_DEBUG == 1)
          {
            cat("CIS_gene_\n")
            cat(str(CIS_gene_INTERVAL))
            cat("\n")
            cat(str(unique(CIS_gene_INTERVAL$ensembl_gene_id)))
            cat("\n")
          }
          
          Block_genes_INTERVAL<-Results_BP_LM_sel_ENSG_sel[is.na(Results_BP_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                             !is.na(Results_BP_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Block_genes_\n")
            cat(str(Block_genes_INTERVAL))
            cat("\n")
            cat(str(unique(Block_genes_INTERVAL$ensembl_gene_id)))
            cat("\n")
          }
          
          #"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS","coefficient_Genotypes_specific_CELL_COUNTS",
          
          df_aggregate<-data.frame()
          
          if(dim(CIS_gene_INTERVAL)[1] >0)
          {
            Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "transcript_id"),
                                                                 which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                 which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                 which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                                 which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"))]
            colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
            colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
            
            if(Condition_DEBUG == 1)
            {
              cat("Block_genes_INTERVAL_subset_0\n")
              cat(str(Block_genes_INTERVAL_subset))
              cat("\n")
              
            }
            
            CIS_gene_INTERVAL_subset<-CIS_gene_INTERVAL[,c(which(colnames(CIS_gene_INTERVAL) == "Cell_Type"),
                                                           which(colnames(CIS_gene_INTERVAL) == "ensembl_gene_id"),
                                                           which(colnames(CIS_gene_INTERVAL) == "HGNC"),
                                                           which(colnames(CIS_gene_INTERVAL) == "CIS_gene_minuslogpvalue"),
                                                           which(colnames(CIS_gene_INTERVAL) == "CIS_Beta"))]
            colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
            colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
            
            if(Condition_DEBUG == 1)
            {
              cat("CIS_gene_INTERVAL_subset_0\n")
              cat(str(CIS_gene_INTERVAL_subset))
              cat("\n")
              
            }
            
            
            Bind_CIS_Block_INTERVAL<-rbind(CIS_gene_INTERVAL_subset,
                                           Block_genes_INTERVAL_subset)
            
            if(Condition_DEBUG == 1)
            {
              cat("Bind_CIS_Block_INTERVAL_0\n")
              cat(str(Bind_CIS_Block_INTERVAL))
              cat("\n")
              
            }
            
            df_aggregate<-rbind(df_aggregate,Bind_CIS_Block_INTERVAL)
            
          }else{
            
            Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "Cell_Type"),
                                                                 which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                 which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                 which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                                 which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"))]
            colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS"
            colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes_specific_CELL_COUNTS"
            
            if(Condition_DEBUG == 1)
            {
              cat("Block_genes_INTERVAL_subset_0\n")
              cat(str(Block_genes_INTERVAL_subset))
              cat("\n")
              
            }
            
            df_aggregate<-rbind(df_aggregate,Block_genes_INTERVAL_subset)
            
          }# dim(CIS_gene_INTERVAL)[1] >0
          
          if(Condition_DEBUG == 1)
          {
            cat("df_aggregate_0\n")
            cat(str(df_aggregate))
            cat("\n")
            
          }
          
         
          
          df_aggregate$VAR<-VAR_sel
          df_aggregate$Significance<-NA
          df_aggregate$RNASeq_source<-NA
          
          if(Condition_DEBUG == 1)
          {
            cat("df_aggregate_1\n")
            cat(str(df_aggregate))
            cat("\n")
            
          }
          
          #### Significance ----
          
          df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3)]<-"YES"
          df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS < 1.3)]<-"NO"
          
          df_aggregate$Significance<-factor(df_aggregate$Significance,
                                            levels=c("NO","YES"),
                                            ordered=T)
          if(Condition_DEBUG == 1)
          {
            cat("df_aggregate$Significance\n")
            cat(sprintf(as.character(names(summary(df_aggregate$Significance)))))
            cat("\n")
            cat(sprintf(as.character(summary(df_aggregate$Significance))))
            cat("\n")
          }
          
          
          
         
          
          df_aggregate_SIG<-df_aggregate[which(df_aggregate$ajusted.minuslogpvalue_Genotypes_specific_CELL_COUNTS >= 1.3),]
          
          if(dim(df_aggregate_SIG)[1] >0)
          {
            if(Condition_DEBUG == 1)
            {
              cat("df_aggregate_SIG_0\n")
              cat(str(df_aggregate_SIG))
              cat("\n")
              
            }
            
            list_per_gene[[z]]<-df_aggregate_SIG
          }#dim(df_aggregate_SIG)[1] >0
        }#z in 1:length(ENSG_array)
        
        if(length(list_per_gene) >0)
        {
          DE_genes_per_VAR = unique(as.data.frame(data.table::rbindlist(list_per_gene, fill=T), stringsAsFactors=F))
          
          if(Condition_DEBUG == 1)
          {
            cat("DE_genes_per_VAR_0\n")
            cat(str(DE_genes_per_VAR))
            cat("\n")
          }
          
          DE_genes_per_VAR.dt<-data.table(DE_genes_per_VAR, key=c("VAR","Cell_Type"))
          
          DE_genes_per_VAR_collapsed<-DE_genes_per_VAR.dt[,.(DE_HGNC_string=paste(HGNC, collapse=";")),by=key(DE_genes_per_VAR.dt)]
          
          if(Condition_DEBUG == 1)
          {
            cat("DE_genes_per_VAR_collapsed_0\n")
            cat(str(DE_genes_per_VAR_collapsed))
            cat("\n")
          }
          
          
          
          DE_genes_per_VAR_collapsed_wide<-as.data.frame(pivot_wider(DE_genes_per_VAR_collapsed,
                                                                     id_cols=c("VAR"),
                                                                     names_from=Cell_Type,
                                                                     values_from=DE_HGNC_string),
                                                         stringsAsFactors=F)
          
          colnames(DE_genes_per_VAR_collapsed_wide)[which(colnames(DE_genes_per_VAR_collapsed_wide) != "VAR")]<-paste(colnames(DE_genes_per_VAR_collapsed_wide)[which(colnames(DE_genes_per_VAR_collapsed_wide) != "VAR")], "DE_HGNC_string",sep="_")
          
          if(Condition_DEBUG == 1)
          {
            cat("DE_genes_per_VAR_collapsed_wide_0:\n")
            cat(str(DE_genes_per_VAR_collapsed_wide))
            cat("\n")
            cat(str(unique(DE_genes_per_VAR_collapsed_wide$VAR)))
            cat("\n")
          }
          # quit(status = 1)
          
          list_DEF[[i]]<-DE_genes_per_VAR_collapsed_wide
          
          # quit(status = 1)
        }#length(list_per_gene) >0
      }#length(ENSG_array) >0
    }#dim(Results_BP_LM_sel)[1] >0
  }#i in 1:length(VARS)
  
  Condition_DEBUG <- 1
  
  if(length(list_DEF) >0)
  {
    
    DTU_gene_df_SIG_collapsed_DEF = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
    
    if(Condition_DEBUG == 1)
    {
      cat("DTU_gene_df_SIG_collapsed_DEF_0\n")
      cat(str(DTU_gene_df_SIG_collapsed_DEF))
      cat("\n")
    }
    
    ZZZ_table<-merge(ZZZ_table,
                     DTU_gene_df_SIG_collapsed_DEF,
                     by="VAR",
                     all.x=T)
    
    if(Condition_DEBUG == 1)
    {
      cat("ZZZ_table_1\n")
      cat(str(ZZZ_table))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Whole_blood_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Whole_blood_DE_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Whole_blood_DTU_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Whole_blood_DTU_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Monocyte_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Monocyte_DE_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Neutrophil_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Neutrophil_DE_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ZZZ_table$Tcell_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ZZZ_table$Tcell_DE_HGNC_string)))))
      cat("\n")
      
      check<-ZZZ_table[(is.na(ZZZ_table$Whole_blood_DE_HGNC_string) &
                          is.na(ZZZ_table$Whole_blood_DTU_HGNC_string) &
                          is.na(ZZZ_table$Monocyte_DE_HGNC_string) &
                          is.na(ZZZ_table$Neutrophil_DE_HGNC_string) &
                          is.na(ZZZ_table$Tcell_DE_HGNC_string)),]
      
      cat("check_1\n")
      cat(str(check))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(check$Whole_blood_DE_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(check$Whole_blood_DE_HGNC_string)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(check$Whole_blood_DTU_HGNC_string))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(check$Whole_blood_DTU_HGNC_string)))))
      cat("\n")
      
      
    }
    
    
    path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
    
    if (file.exists(path7)){
      
      
      
      
    } else {
      dir.create(file.path(path7))
      
    }
    
    if(Condition_DEBUG == 1)
    {
      cat("path7:\t")
      cat(sprintf(as.character(path7)))
      cat("\n")
    }
    
    setwd(path7)
    
    write.table(ZZZ_table, sep="\t", quote=F, row.names = F,file="ZZZ_in_process.tsv")
    
  }# length(list_DEF) >0
  
}

ADD_Proxys_CSQ = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  
  #### CURATED TABLE ----
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  setwd(path7)
  
  ZZZ_table<-as.data.frame(fread(file="ZZZ_in_process.tsv", sep ="\t", header =T), stringsAsFactors=F)
  
  
  
  cat("ZZZ_table_0:\n")
  cat(str(ZZZ_table))
  cat("\n")
  cat(str(unique(ZZZ_table$VAR)))
  cat("\n")
  
  ### Read Proxy_CSQ_TABLE----
  
  
  
  Proxy_CSQ_TABLE<-unique(as.data.frame(fread(file=opt$Proxy_CSQ_TABLE, sep="\t", header=T) , stringsAsFactors=T))
  
  
  cat("Proxy_CSQ_TABLE\n")
  cat(str(Proxy_CSQ_TABLE))
  cat("\n")
  cat(str(unique(Proxy_CSQ_TABLE$VAR)))
  cat("\n")
  
  ### Read Proxy_R2_TABLE----
  
  
  
  Proxy_R2_TABLE<-unique(as.data.frame(fread(file=opt$Proxy_R2_TABLE, sep="\t", header=T) , stringsAsFactors=T))
  
  Proxy_R2_TABLE$R2<-round(Proxy_R2_TABLE$R2,2)
  
  cat("Proxy_R2_TABLE\n")
  cat(str(Proxy_R2_TABLE))
  cat("\n")
  cat(str(unique(Proxy_R2_TABLE$VAR)))
  cat("\n")
  
  #### Merge Proxys----
  
  Proxy_merge<-merge(Proxy_R2_TABLE,
                     Proxy_CSQ_TABLE,
                     by=colnames(Proxy_R2_TABLE)[which(colnames(Proxy_R2_TABLE)%in%colnames(Proxy_CSQ_TABLE))],
                     all=T)
  
  cat("Proxy_merge\n")
  cat(str(Proxy_merge))
  cat("\n")
  cat(str(unique(Proxy_merge$VAR)))
  cat("\n")
  cat(str(unique(Proxy_merge$Proxy_VAR)))
  cat("\n")
  
  #### data.table to string ----
  Condition_DEBUG <- 0
  
  Proxy_merge.dt<-data.table(Proxy_merge, key=c("VAR","Proxy_VAR","rsid","Proxy_rsid","R2","HGNC"))
  
  Collapse_CSQ<-as.data.frame(Proxy_merge.dt[,.(CSQ_string=paste(unique(VEP_DEF_LABELS), collapse=",")), by=key(Proxy_merge.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Collapse_CSQ\n")
    cat(str(Collapse_CSQ))
    cat("\n")
    cat(str(unique(Collapse_CSQ$VAR)))
    cat("\n")
    cat(str(unique(Collapse_CSQ$Proxy_VAR)))
    cat("\n")
  }
  Collapse_CSQ.dt<-data.table(Collapse_CSQ, key=c("VAR","Proxy_VAR","rsid","Proxy_rsid","R2"))
  
  Collapse_HGNC<-as.data.frame(Collapse_CSQ.dt[,.(HGNC_string=paste(paste(HGNC,'(',CSQ_string,')',sep=''), collapse=';')), by=key(Collapse_CSQ.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Collapse_HGNC\n")
    cat(str(Collapse_HGNC))
    cat("\n")
    cat(str(unique(Collapse_HGNC$VAR)))
    cat("\n")
    cat(str(unique(Collapse_HGNC$Proxy_VAR)))
    cat("\n")
  }
  
  Collapse_HGNC.dt<-data.table(Collapse_HGNC, key=c("VAR"))
  
  Collapse_Proxy_rsid<-as.data.frame(Collapse_HGNC.dt[,.(Proxy_rsid_string=paste(paste(Proxy_rsid,
                                                                                       R2,
                                                                                       HGNC_string,
                                                                                       sep="__"), collapse = "|")), by=key(Collapse_HGNC.dt)], stringsAsFactors=F)
  if(Condition_DEBUG == 1)
  {
    cat("Collapse_Proxy_rsid\n")
    cat(str(Collapse_Proxy_rsid))
    cat("\n")
    cat(str(unique(Collapse_Proxy_rsid$VAR)))
    cat("\n")
  }
  
  ##### LAST MERGE ----
  
  ZZZ_table<-merge(ZZZ_table,
                   Collapse_Proxy_rsid,
                   by="VAR",
                   all.x=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_1\n")
    cat(str(ZZZ_table))
    cat("\n")
  }
  
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  if(Condition_DEBUG == 1)
  {
    cat("path7:\t")
    cat(sprintf(as.character(path7)))
    cat("\n")
  }
  
  setwd(path7)
  
  write.table(ZZZ_table, sep="\t", quote=F, row.names = F,file="ZZZ_in_process.tsv")
  
}

ADD_Screenings = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  
  #### CURATED TABLE ----
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  setwd(path7)
  
  ZZZ_table<-as.data.frame(fread(file="ZZZ_in_process.tsv", sep ="\t", header =T), stringsAsFactors=F)
  
  
  
  cat("ZZZ_table_0:\n")
  cat(str(ZZZ_table))
  cat("\n")
  cat(str(unique(ZZZ_table$VAR)))
  cat("\n")
  
  Condition_DEBUG <-0
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  
  cat("CUMMULATIVE_CLASSES_0:\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
  cat("\n")
  
  CUMMULATIVE_CLASSES$comparison_VAR<-gsub("^chr","",CUMMULATIVE_CLASSES$VAR)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES\n")
    cat(str(CUMMULATIVE_CLASSES))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted<-CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$carried_variants == CUMMULATIVE_CLASSES$comparison_VAR),]
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted<-droplevels(CUMMULATIVE_CLASSES_restricted[which(CUMMULATIVE_CLASSES_restricted$Cell_Type == "ALL_CT"),])
  
  CUMMULATIVE_CLASSES_restricted$enhancer_classif<-"NA"
  CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif<-"NA"
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_1\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted$enhancer_classif[which(CUMMULATIVE_CLASSES_restricted$enhancer_CLASS_TILES == 0)]<-"NO_TILES_with_enhancer_activity"
  CUMMULATIVE_CLASSES_restricted$enhancer_classif[which(CUMMULATIVE_CLASSES_restricted$enhancer_CLASS_TILES > 0)]<-"AT_LEAST_1_TILE_with_enhancer_activity"
  
  CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif[which(CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_CLASS_TILES == 0)]<-"NO_TILES_with_E_Plus_ASE_activity"
  CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif[which(CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_CLASS_TILES > 0)]<-"AT_LEAST_1_TILE_with_E_Plus_ASE_activity"
  
  CUMMULATIVE_CLASSES_restricted$interaction_5<-interaction(CUMMULATIVE_CLASSES_restricted$enhancer_classif,
                                                            CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif,
                                                            sep="|",lex.order = T)
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_2\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat("E_Plus_ASE_CLASS_TILES\n")
    cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_CLASS_TILES)))))
    cat("\n")
    cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_CLASS_TILES))))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR))) # 94
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(CUMMULATIVE_CLASSES_restricted$enhancer_classif))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(CUMMULATIVE_CLASSES_restricted$enhancer_classif)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(CUMMULATIVE_CLASSES_restricted$interaction_5))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(CUMMULATIVE_CLASSES_restricted$interaction_5)))))
    cat("\n")
  }
  
  
  CUMMULATIVE_CLASSES_restricted$MPRA_CLASS<-"NA"
  
  CUMMULATIVE_CLASSES_restricted$MPRA_CLASS[which(CUMMULATIVE_CLASSES_restricted$interaction_5 == "AT_LEAST_1_TILE_with_enhancer_activity|AT_LEAST_1_TILE_with_E_Plus_ASE_activity")]<-"AT_LEAST_1_TILE_with_E_Plus_ASE_activity"
  CUMMULATIVE_CLASSES_restricted$MPRA_CLASS[which(CUMMULATIVE_CLASSES_restricted$interaction_5 == "AT_LEAST_1_TILE_with_enhancer_activity|NO_TILES_with_E_Plus_ASE_activity")]<-"AT_LEAST_1_TILE_with_enhancer_activity"
  CUMMULATIVE_CLASSES_restricted$MPRA_CLASS[which(CUMMULATIVE_CLASSES_restricted$interaction_5 == "NO_TILES_with_enhancer_activity|AT_LEAST_1_TILE_with_E_Plus_ASE_activity")]<-"AT_LEAST_1_TILE_with_E_Plus_ASE_activity" # doesn't exist
  CUMMULATIVE_CLASSES_restricted$MPRA_CLASS[which(CUMMULATIVE_CLASSES_restricted$interaction_5 == "NO_TILES_with_enhancer_activity|NO_TILES_with_E_Plus_ASE_activity")]<-"NO_enhancer_activity"
  
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_3\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
  }
  
  
  indx.int<-c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
              which(colnames(CUMMULATIVE_CLASSES_restricted) == "MPRA_CLASS"))
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[,indx.int])
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
  }
  
  ### genIE_input regulation #----
  
  
  genIE_input<-as.data.frame(fread(file=opt$genIE_input, sep =",", header =T), stringsAsFactors=F)
  
  
 
  cat("genIE_input_0:\n")
  cat(str(genIE_input))
  cat("\n")
  cat(str(unique(genIE_input$VAR)))
  cat("\n")
  
  genIE_input_subset<-unique(genIE_input[which(genIE_input$HGNC != ""),])
  
  if(Condition_DEBUG == 1)
  {
    cat("genIE_input_subset_0\n")
    cat(str(genIE_input_subset))
    cat("\n")
    cat(sprintf(as.character(genIE_input_subset$HGNC)))
    cat("\n")
  }
  
  
  HGNC_tested<-c("FOXP1","CUX1","TNRC6A","SH2B3","BRAP","UGCG","BID","NBN","C2CD5","EPB41")
  
  genIE_input_subset<-unique(genIE_input_subset[which(genIE_input_subset$HGNC%in%HGNC_tested),])
  
  if(Condition_DEBUG == 1)
  {
    cat("genIE_input_subset_1\n")
    cat(str(genIE_input_subset))
    cat("\n")
    cat(sprintf(as.character(genIE_input_subset$HGNC)))
    cat("\n")
  }
  
  genIE_input_subset<-unique(genIE_input_subset[,c(which(colnames(genIE_input_subset) == "VAR"),
                                                   which(colnames(genIE_input_subset) == "HGNC"))])
  if(Condition_DEBUG == 1)
  {
    cat("genIE_input_subset_2\n")
    cat(str(genIE_input_subset))
    cat("\n")
  }
  
  ### genIE_active regulation #----
  
  
  genIE_active<-as.data.frame(fread(file=opt$genIE_active, sep ="\t", header =T), stringsAsFactors=F)
  
  
  
  cat("genIE_active_0:\n")
  cat(str(genIE_active))
  cat("\n")
  cat(str(unique(genIE_active$VAR)))
  cat("\n")
  
  #### genIE CLASSIFICATION ----
  
  genIE_input_subset$genIE_CLASS<-"NA"
  
  
  genIE_input_subset$genIE_CLASS[which(genIE_input_subset$VAR%in%genIE_active$VAR)]<-"genIE_ACTIVE"
  genIE_input_subset$genIE_CLASS[-which(genIE_input_subset$VAR%in%genIE_active$VAR)]<-"genIE_INACTIVE"
  
  
  
  
  genIE_input_subset<-genIE_input_subset[,-which(colnames(genIE_input_subset) == "HGNC")]
  
  if(Condition_DEBUG == 1)
  {
    cat("genIE_input_subset_3\n")
    cat(str(genIE_input_subset))
    cat("\n")
    cat(str(unique(genIE_input_subset$VAR)))
    cat("\n")
  }
  ##### Manual curation QTLS ----
  
  VARS_OT_QTLS<-c("chr4_1008212_C_T","chr9_135857646_C_T","chr7_50444152_G_T","chr1_202129205_G_A","chr15_65174438_C_A","chr13_28604007_T_C","chr6_34947254_A_G","chr9_114663385_T_C","chr17_56603493_C_T","chr6_41924998_C_T","chr6_41952511_T_G","chr5_35476470_G_T","chr22_50949811_T_C","chr7_100309180_A_G","chr1_92925654_G_C","chr1_29217311_G_A","chr2_219020958_C_T","chr16_155132_C_T","chr7_99760955_G_C","chr1_198680015_G_A","chr12_111844956_C_T","chr15_64349614_G_A")
  CLASSIF<-c("eQTL_OT_ABSENT","eQTL_OT_ABSENT","eQTL_OT_ABSENT","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_ABSENT","eQTL_OT_FOUND","eQTL_OT_ABSENT","eQTL_OT_FOUND","eQTL_OT_ABSENT","eQTL_OT_ABSENT","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_NOT_TESTED","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_FOUND","eQTL_OT_ABSENT")
  
  df<-as.data.frame(cbind(VARS_OT_QTLS,CLASSIF), stringsAsFactors=F)
  
  colnames(df)<-c("VAR","Replication_OT_QTL")
  
  ##### LAST MERGE ----
  
  ZZZ_table<-merge(ZZZ_table,
                   CUMMULATIVE_CLASSES_restricted_subset,
                   by="VAR",
                   all.x=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_1\n")
    cat(str(ZZZ_table))
    cat("\n")
  }
  
  ZZZ_table<-merge(ZZZ_table,
                   genIE_input_subset,
                   by="VAR",
                   all.x=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_2\n")
    cat(str(ZZZ_table))
    cat("\n")
  }
  
  ZZZ_table<-merge(ZZZ_table,
                   df,
                   by="VAR",
                   all.x=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_3\n")
    cat(str(ZZZ_table))
    cat("\n")
  }
  
  ZZZ_table$MPRA_CLASS[is.na(ZZZ_table$MPRA_CLASS)]<-"NOT_SCREENED_MPRA"
  ZZZ_table$genIE_CLASS[is.na(ZZZ_table$genIE_CLASS)]<-"NOT_SCREENED_genIE"
  
  indx<-which(ZZZ_table$Mechanistic_Class == "No_RNA_Seq_HET_carriers")
  
  cat("indx:\n")
  cat(str(indx))
  cat("\n")
  
  ZZZ_table$Replication_OT_QTL[indx]<-"No_RNA_Seq_HET_carriers"
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_4\n")
    cat(str(ZZZ_table))
    cat("\n")
  }
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  if(Condition_DEBUG == 1)
  {
    cat("path7:\t")
    cat(sprintf(as.character(path7)))
    cat("\n")
  }
  
  setwd(path7)
  
  write.table(ZZZ_table, sep="\t", quote=F, row.names = F,file="ZZZ_in_process.tsv")
  
}


ADD_phenotypes = function(option_list)
{
  
  suppressMessages(library("BiocGenerics", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("S4Vectors", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("IRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("Biobase", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GenomicFeatures", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("OrganismDbi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("GO.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("org.Hs.eg.db", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("Homo.sapiens", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("gwascat", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("rtracklayer", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  suppressMessages(library("liftOver",lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  
  #### CURATED TABLE ----
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  setwd(path7)
  
  ZZZ_table<-as.data.frame(fread(file="ZZZ_in_process.tsv", sep ="\t", header =T), stringsAsFactors=F)
  
  
  
  cat("ZZZ_table_0:\n")
  cat(str(ZZZ_table))
  cat("\n")
  cat(str(unique(ZZZ_table$VAR)))
  cat("\n")
  
  Condition_DEBUG <-0
  
  #### READ and transform Prob_Threshold ----
  
  Prob_Threshold = opt$Prob_Threshold
  
  cat("Prob_Threshold\n")
  cat(sprintf(as.character(Prob_Threshold)))
  cat("\n")
  
  #### CURATED TABLE ----
  
  AF_and_CSQ_file<-readRDS(file=opt$AF_and_CSQ_file)
  
  
  cat("AF_and_CSQ_file_0:\n")
  cat(str(AF_and_CSQ_file))
  cat("\n")
  
  indx.int<-c(which(colnames(AF_and_CSQ_file) == "VAR"),which(colnames(AF_and_CSQ_file) == "rs"),which(colnames(AF_and_CSQ_file) == "maf_origin"),which(colnames(AF_and_CSQ_file) == "VEP_DEF_LABELS_wCSQ"))
  
  AF_and_CSQ_file_subset<-AF_and_CSQ_file[,indx.int]
  
  if(Condition_DEBUG == 1)
  {
    cat("AF_and_CSQ_file_subset_0:\n")
    cat(str(AF_and_CSQ_file_subset))
    cat("\n")
  }
  
  
  
  
  #### ALL_dB ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB, sep ="\t", header =T), stringsAsFactors=F)
  
  ALL_dB$Lineage<-"NA"
  
  cat("ALL_dB_0:\n")
  cat(str(ALL_dB))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ALL_dB$phenotype))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ALL_dB$phenotype)))))
  cat("\n")
  
  
  ##### P vs NP classification -----
  
  # lymph lymph_p wbc
  
  
  erythroid.lineage <- c("mchc","ret_p","ret","mrv","mscv","mch","mcv","rdw_cv","MicroR","hgb","rbc","hct","irf","hlr_p","hlr")
  
  
  ALL_dB$Lineage[which(ALL_dB$phenotype%in%erythroid.lineage)]<-"erythroid_lineage"
  
  megakaryocytic.lineage <- c("plt","pct","mpv","pdw")
  
  
  ALL_dB$Lineage[which(ALL_dB$phenotype%in%megakaryocytic.lineage)]<-"megakaryocytic_lineage"
  
  granulocyte_monocyte.lineage <- c("neut","mono","baso","eo") 
  
  
  ALL_dB$Lineage[which(ALL_dB$phenotype%in%granulocyte_monocyte.lineage)]<-"granulocyte_monocyte_lineage"
  
  lymphocyte.lineage <- c("lymph") 
  
  ALL_dB$Lineage[which(ALL_dB$phenotype%in%lymphocyte.lineage)]<-"lymphocyte_lineage"
  
  # wbc.lineage <- c()
  # 
  # 
  # ALL_dB$Lineage[which(ALL_dB$phenotype%in%wbc.lineage)]<-"wbc_mix"
  
  # "wbc_mix"
  
  
  
  ALL_dB$Lineage<-factor(as.character(ALL_dB$Lineage),
                         levels=c("erythroid_lineage","megakaryocytic_lineage","granulocyte_monocyte_lineage","lymphocyte_lineage"), ordered = T)
  
  cat(sprintf(as.character(names(summary(ALL_dB$Lineage)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_dB$Lineage))))
  cat("\n")
  
  
  
  #### Filter for Prob threshold ----
  
  ALL_dB_filter_Prob<-ALL_dB[which(ALL_dB$finemap_prob >= Prob_Threshold),]
  
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_\n")
    str(ALL_dB_filter_Prob)
    cat("\n")
    cat(str(unique(ALL_dB_filter_Prob$VAR)))
    cat("\n") #
  }
  
  #### Find which ons are in the Supp 4 table ----
  
  
  ALL_dB_filter_Prob_subset<-ALL_dB_filter_Prob[which(ALL_dB_filter_Prob$VAR%in%ZZZ_table$VAR),]
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_0\n")
    str(ALL_dB_filter_Prob_subset)
    cat("\n")
    cat(str(unique(ALL_dB_filter_Prob_subset$VAR)))
    cat("\n") #
    cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_Prob_subset$Lineage))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_Prob_subset$Lineage)))))
    cat("\n")
  }
  
  ALL_dB_filter_Prob_subset$phenotype_DEF<-revalue(ALL_dB_filter_Prob_subset$phenotype, 
                                                   c("plt"="PLT#",
                                                     "mpv"="MPV",
                                                     "pdw"="PDW",
                                                     "pct"="PCT",
                                                     "H-IPF"="H-IPF",
                                                     "IPF_perc"="IPF%",
                                                     "IPF"="IPF#",
                                                     "P-LCR"="P-LCR",
                                                     "rbc"="RBC#",
                                                     "mcv"="MCV",
                                                     "hct"="HCT",
                                                     "mch"="MCH",
                                                     "mchc"="MCHC",
                                                     "hgb"="HGB",
                                                     "rdw_cv"='RDW',
                                                     "MacrorR"='MacroR',
                                                     "MicroR"='MicroR',
                                                     "RBC-He"='RBC-He',
                                                     "Delta-He"='Delta-He',
                                                     "Hyer-He"='Hyper-He',
                                                     "RDW-SD"='RDW-SD',
                                                     "RPI"='RPI',
                                                     "HFR"='HFR',
                                                     "MFR"='MFR',
                                                     "LFR"='LFR',
                                                     "IG"='IG#',
                                                     "IG_perc"='IG%',
                                                     "ret"='RET#',
                                                     "ret_p"='RET%',
                                                     "irf"='IRF',
                                                     "hlr"='HLSR#',
                                                     "hlr_p"='HLSR%',
                                                     "mrv"='MRV',
                                                     "mscv"='MSCV',
                                                     "RET-FSC"='RET-FSC',
                                                     "RET-RBC-FSC"='RET-FSC',
                                                     "IRF-FSC"='IRF-FSC',
                                                     "RET-He"='RET-He',
                                                     "RET-UPP"='RET-UPP',
                                                     "lymph"='LYMPH#',
                                                     "lymph_p"='LYMPH%',
                                                     "LY-FSC"="LY-FSC",
                                                     "LY-FSC-DW"="LY-FSC-DW",
                                                     "LY-SSC"="LY-SSC",
                                                     "LY-SSC-DW"="LY-SSC-DW",
                                                     "LY-SFL"="LY-SFL",
                                                     "LY-SFL-DW"="LY-SFL-DW",
                                                     "baso"='BASO#',
                                                     "baso_p"='BASO%',
                                                     "eo"='EO#',
                                                     "eo_p"='EO%',
                                                     "mono"='MONO#',
                                                     "mono_p"='MONO%',
                                                     "MO-FSC"="MO-FSC",
                                                     "MO-FSC-DW"="MO-FSC-DW",
                                                     "MO-SSC"="MO-SSC",
                                                     "MO-SSC-DW"="MO-SSC-DW",
                                                     "MO-SFL"="MO-SFL",
                                                     "MO-SFL-DW"="MO-SFL-DW",
                                                     "neut"='NEUT#',
                                                     "neut_p"='NEUT%',
                                                     "NE-FSC"="NE-FSC",
                                                     "NE-FSC-DW"="NE-FSC-DW",
                                                     "NE-SSC"="NE-SSC",
                                                     "NE-SSC-DW"="NE-SSC-DW",
                                                     "NE-SFL"="NE-SFL",
                                                     "NE-SFL-DW"="NE-SFL-DW",
                                                     "wbc"='WBC#'))
  
  ALL_dB_filter_Prob_subset$phenotype_DEF<-factor(as.character(ALL_dB_filter_Prob_subset$phenotype_DEF),
                                                  levels=c("RBC#","MCV","HCT","MCH","MCHC","HGB",'RDW','MacroR','MicroR','RBC-He','Delta-He','Hyper-He','RDW-SD','RPI','HFR','MFR','LFR','IG#','IG%',
                                                           'RET#','RET%','IRF','HLSR#','HLSR%','MRV','MSCV',"RET-FSC","RET-RBC-FSC","IRF-FSC","RET-He","RET-UPP",
                                                           "PLT#","MPV","PDW","PCT","H-IPF","IPF%","IPF#","P-LCR",
                                                           'MONO#','MONO%',"MO-FSC","MO-FSC-DW","MO-SSC","MO-SSC-DW","MO-SFL","MO-SFL-DW",
                                                           'NEUT#','NEUT%',"NE-FSC","NE-FSC-DW","NE-SSC","NE-SSC-DW","NE-SFL","NE-SFL-DW",
                                                           'BASO#','BASO%','EO#','EO%',
                                                           'LYMPH#','LYMPH%',"LY-FSC","LY-FSC-DW","LY-SSC","LY-SSC-DW","LY-SFL","LY-SFL-DW",
                                                           'WBC#'), ordered = T)
  
  # quit(status = 1)
  
  
  
  ALL_dB_filter_Prob_subset<-droplevels(ALL_dB_filter_Prob_subset[order(ALL_dB_filter_Prob_subset$phenotype_DEF),])
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_1\n")
    str(ALL_dB_filter_Prob_subset)
    cat("\n")
    cat(str(unique(ALL_dB_filter_Prob_subset$VAR)))
    cat("\n") #
  }
  
  ALL_dB_filter_Prob_subset.dt<-data.table(ALL_dB_filter_Prob_subset, key=c("rs","VAR"))
  
  
  # cat("ALL_dB_filter_Prob_subset.dt\n")
  # cat(str(ALL_dB_filter_Prob_subset.dt))
  # cat("\n")
  
  ALL_dB_filter_Prob_subset_compact<-as.data.frame(ALL_dB_filter_Prob_subset.dt[,.(phenotype_string=paste(phenotype, collapse=";"),
                                                                                   phenotype_DEF_string=paste(phenotype, collapse=";"),
                                                                                   Lineage_string=paste(unique(Lineage[!is.na(Lineage)]), collapse=";")), 
                                                                                by=key(ALL_dB_filter_Prob_subset.dt)], stringsAsFactors=F)
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_compact_0\n")
    cat(str(ALL_dB_filter_Prob_subset_compact))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_Prob_subset_compact$Lineage_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_Prob_subset_compact$Lineage_string)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_Prob_subset_compact$phenotype_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_Prob_subset_compact$phenotype_string)))))
    cat("\n")
  }
  
  
  #### RESCUE MISS ASSIGNED----
  
  
  correlated_phenotypes<-c("mono_p","neut_p","baso_p","eo_p","lymph_p","wbc")
  
  indx.correlated<-grep(paste(correlated_phenotypes,collapse ="|"), ALL_dB_filter_Prob_subset_compact$phenotype_string)
  
  Correlated_phenotypes<-ALL_dB_filter_Prob_subset_compact[indx.correlated,]
  
  if(Condition_DEBUG == 1)
  {
    cat("Correlated_phenotypes\n")
    cat(str(Correlated_phenotypes))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Correlated_phenotypes$Lineage_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Correlated_phenotypes$Lineage_string)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Correlated_phenotypes$phenotype_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Correlated_phenotypes$phenotype_string)))))
    cat("\n")
  }
  
  
  MISS_asigned_categories<-c("baso_p;mscv;pct;plt","mono_p")
  
  MISS_ASIGNED<-ALL_dB_filter_Prob_subset_compact[which(ALL_dB_filter_Prob_subset_compact$phenotype_string %in% MISS_asigned_categories),]
  
  if(Condition_DEBUG == 1)
  {
    cat("MISS_ASIGNED_0\n")
    cat(str(MISS_ASIGNED))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(MISS_ASIGNED$Lineage_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(MISS_ASIGNED$Lineage_string)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(MISS_ASIGNED$phenotype_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(MISS_ASIGNED$phenotype_string)))))
    cat("\n")
  }
  
  MISS_ASIGNED$Lineage_string[(MISS_ASIGNED$phenotype_string == "baso_p;mscv;pct;plt")]<-"erythroid_lineage;megakaryocytic_lineage;granulocyte_monocyte_lineage"
  
  MISS_ASIGNED$Lineage_string[(MISS_ASIGNED$phenotype_string == "mono_p")]<-"granulocyte_monocyte_lineage"
  
  if(Condition_DEBUG == 1)
  {
    cat("MISS_ASIGNED_1\n")
    cat(str(MISS_ASIGNED))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(MISS_ASIGNED$Lineage_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(MISS_ASIGNED$Lineage_string)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(MISS_ASIGNED$phenotype_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(MISS_ASIGNED$phenotype_string)))))
    cat("\n")
  }
  
  ALL_dB_filter_Prob_subset_compact<-ALL_dB_filter_Prob_subset_compact[-which(ALL_dB_filter_Prob_subset_compact$VAR%in%MISS_ASIGNED$VAR),]
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_compact_1\n")
    cat(str(ALL_dB_filter_Prob_subset_compact))
    cat("\n")
  }
  
  ALL_dB_filter_Prob_subset_compact<-rbind(MISS_ASIGNED,ALL_dB_filter_Prob_subset_compact)
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_compact_2\n")
    cat(str(ALL_dB_filter_Prob_subset_compact))
    cat("\n")
  }
  
  #### grep multiple categories ----
  
  indx_Multi_Lineage<-grep(";",ALL_dB_filter_Prob_subset_compact$Lineage_string)
  
  
  Multi_Lineage_df<-ALL_dB_filter_Prob_subset_compact[indx_Multi_Lineage,]
  
  if(Condition_DEBUG == 1)
  {
    cat("Multi_Lineage_df_1\n")
    cat(str(Multi_Lineage_df))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Multi_Lineage_df$Lineage_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Multi_Lineage_df$Lineage_string)))))
    cat("\n")
  }
  
  Lineage_restricted_df<-ALL_dB_filter_Prob_subset_compact[-indx_Multi_Lineage,]
  
  if(Condition_DEBUG == 1)
  {
    cat("Lineage_restricted_df_1\n")
    cat(str(Lineage_restricted_df))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Lineage_restricted_df$Lineage_string))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Lineage_restricted_df$Lineage_string)))))
    cat("\n")
  }
  
  
  
  ALL_dB_filter_Prob_subset_compact$Multi_Lineage<-"NA"
  
  ALL_dB_filter_Prob_subset_compact$Multi_Lineage[which(ALL_dB_filter_Prob_subset_compact$VAR%in%Multi_Lineage_df$VAR)]<-"Multi_Lineage"
  
  ALL_dB_filter_Prob_subset_compact$Multi_Lineage[which(ALL_dB_filter_Prob_subset_compact$VAR%in%Lineage_restricted_df$VAR)]<-"Lineage_restricted"
  
  ALL_dB_filter_Prob_subset_compact$Multi_Lineage<-factor(as.character(ALL_dB_filter_Prob_subset_compact$Multi_Lineage),
                                                          levels=c("Lineage_restricted","Multi_Lineage"), ordered = T)
  
  ALL_dB_filter_Prob_subset_compact<-ALL_dB_filter_Prob_subset_compact[,-which(colnames(ALL_dB_filter_Prob_subset_compact) == "phenotype_string")]
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_compact_3\n")
    cat(str(ALL_dB_filter_Prob_subset_compact))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_Prob_subset_compact$Multi_Lineage))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_Prob_subset_compact$Multi_Lineage)))))
    cat("\n")
  }
  
  
  check<-ALL_dB_filter_Prob_subset_compact[is.na(ALL_dB_filter_Prob_subset_compact$Multi_Lineage),]
  
  if(Condition_DEBUG == 1)
  {
    cat("check_1\n")
    cat(str(check))
    cat("\n")
  }
  # cat(sprintf(as.character(names(summary(as.factor(check$phenotype_DEF_string))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(check$phenotype_DEF_string)))))
  # cat("\n")
  
  
  # quit(status = 1)
  
  
  
  
  ALL_dB_filter_Prob_subset_compact$chr<-gsub("_.+$","",ALL_dB_filter_Prob_subset_compact$VAR)
  ALL_dB_filter_Prob_subset_compact$pos37<-gsub("^chr[^_]+_","",ALL_dB_filter_Prob_subset_compact$VAR)
  ALL_dB_filter_Prob_subset_compact$pos37<-gsub("_.+$","",ALL_dB_filter_Prob_subset_compact$pos37)
  ALL_dB_filter_Prob_subset_compact$ref<-gsub("^chr[^_]+_[^_]+_","",ALL_dB_filter_Prob_subset_compact$VAR)
  ALL_dB_filter_Prob_subset_compact$ref<-gsub("_.+$","",ALL_dB_filter_Prob_subset_compact$ref)
  ALL_dB_filter_Prob_subset_compact$alt<-gsub("^chr[^_]+_[^_]+_[^_]+_","",ALL_dB_filter_Prob_subset_compact$VAR)
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_compact_2\n")
    cat(str(ALL_dB_filter_Prob_subset_compact))
    cat("\n")
  }
  
  # cat(sprintf(as.character(names(summary(as.factor(ALL_dB_filter_Prob_subset$phenotype_DEF))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(ALL_dB_filter_Prob_subset$phenotype_DEF)))))
  # cat("\n")
  
  
  
  ###### LiftOver 37 -> 38 ----
  
  
  gr_VARS <- GRanges(
    seqnames = as.character(gsub("chr","",ALL_dB_filter_Prob_subset_compact$chr)),
    ranges=IRanges(
      start=as.numeric(ALL_dB_filter_Prob_subset_compact$pos37),
      end=as.numeric(ALL_dB_filter_Prob_subset_compact$pos37),
      name=ALL_dB_filter_Prob_subset_compact$VAR))
  
  # cat("gr_VARS\n")
  # str(gr_VARS)
  # cat("\n")
  
  VAR_df<-data.frame(chr=as.character(paste('chr',seqnames(gr_VARS), sep='')),
                     pos37=start(gr_VARS),
                     ref=ALL_dB_filter_Prob_subset_compact$ref,
                     alt=ALL_dB_filter_Prob_subset_compact$alt,
                     VAR=ALL_dB_filter_Prob_subset_compact$VAR,
                     stringsAsFactors = F)
  
  if(Condition_DEBUG == 1)
  {
    cat("VAR_df_\n")
    str(VAR_df)
    cat("\n")
  }
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/nfs/team151/software/manuel_R_ext_data_4_1/hg19ToHg38.over.chain")
  
  seqlevelsStyle(gr_VARS) = "UCSC"  # necessary
  gr_VARS38 = liftOver(gr_VARS, ch)
  gr_VARS38 = unlist(gr_VARS38)
  genome(gr_VARS38) = "hg38"
  
  if(length(gr_VARS38) >0)
  {
    
    chr_38<-as.character(seqnames(gr_VARS38))
    names_38<-as.character(names(gr_VARS38))
    
    ref_VAR38<-gsub("^chr[^_]+_[0-9]+_","",names_38)
    ref_VAR38<-gsub("_.+$","",ref_VAR38)
    
    
    # cat("ref_VAR38\n")
    # cat(sprintf(as.character(ref_VAR38)))
    # cat("\n")
    
    alt_VAR38<-gsub("^chr[^_]+_[0-9]+_[^_]+_","",names_38)
    # alt_VAR38<-gsub("_.+$","",alt_VAR38)
    
    
    # cat("alt_VAR38\n")
    # cat(sprintf(as.character(alt_VAR38)))
    # cat("\n")
    
    
    
    
    VAR_38_df<-data.frame(chr=as.character(seqnames(gr_VARS38)),
                          pos_38=start(gr_VARS38),
                          ref=ref_VAR38,
                          alt=alt_VAR38,
                          VAR=names(gr_VARS38),
                          stringsAsFactors = F)
    
    VAR_38_df$VAR_38<-paste(VAR_38_df$chr,VAR_38_df$pos_38,VAR_38_df$ref,VAR_38_df$alt,sep='_')
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_38_df_1\n")
      str(VAR_38_df)
      cat("\n")
    }
    
    
    VAR_ALL_dB_filter_Prob_subset_compact_df<-unique(merge(VAR_df,
                                                           VAR_38_df,
                                                           by=c("chr","ref","alt","VAR"),
                                                           all=T))
    
    VAR_ALL_dB_filter_Prob_subset_compact_df$VAR_38[is.na(VAR_ALL_dB_filter_Prob_subset_compact_df$VAR_38)]<-"ABSENT"
    
    if(Condition_DEBUG == 1)
    {
      cat("VAR_ALL_dB_filter_Prob_subset_compact_df_2\n")
      str(VAR_ALL_dB_filter_Prob_subset_compact_df)
      cat("\n")
    }
    
    check.ABSENT<-VAR_ALL_dB_filter_Prob_subset_compact_df[which(VAR_ALL_dB_filter_Prob_subset_compact_df$VAR_38 == "ABSENT"),]
    #
    if(Condition_DEBUG == 1)
    {
      cat("check.ABSENT\n")
      str(check.ABSENT) #7
      cat("\n")
    }
    
    
    
    ALL_dB_filter_Prob_subset_compact<-unique(merge(ALL_dB_filter_Prob_subset_compact,
                                                    VAR_ALL_dB_filter_Prob_subset_compact_df[which(VAR_ALL_dB_filter_Prob_subset_compact_df$VAR_38 != "ABSENT"),],
                                                    by=c("VAR","chr","pos37","ref","alt")))
    if(Condition_DEBUG == 1)
    {
      cat("df_2\n")
      str(df)
      cat("\n")
    }
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_VARS38) >0
  
  
  
  # ALL_dB_filter_Prob_subset_compact$chr.pos_ref_alt<-gsub("^chr","",ALL_dB_filter_Prob_subset_compact$VAR_38)
  # ALL_dB_filter_Prob_subset_compact$chr.pos_ref_alt<-sub("_",":",ALL_dB_filter_Prob_subset_compact$chr.pos_ref_alt)
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_filter_Prob_subset_compact\n")
    cat(str(ALL_dB_filter_Prob_subset_compact))
    cat("\n")
    cat(str(unique(ALL_dB_filter_Prob_subset_compact$VAR)))
    cat("\n")
  }
  
  
  indx.select<-c(which(colnames(ALL_dB_filter_Prob_subset_compact) == "VAR"),which(colnames(ALL_dB_filter_Prob_subset_compact) == "rs"),which(colnames(ALL_dB_filter_Prob_subset_compact) == "VAR_38"),
                 which(colnames(ALL_dB_filter_Prob_subset_compact) == "Multi_Lineage"),which(colnames(ALL_dB_filter_Prob_subset_compact) == "Lineage_string"),which(colnames(ALL_dB_filter_Prob_subset_compact) == "phenotype_DEF_string"))
  
  
  
  ALL_dB_FINAL<-unique(ALL_dB_filter_Prob_subset_compact[,indx.select])
  
  if(Condition_DEBUG == 1)
  {
    cat("ALL_dB_FINAL\n")
    cat(str(ALL_dB_FINAL))
    cat("\n")
    cat(str(unique(ALL_dB_FINAL$VAR)))
    cat("\n")
  }
  
  
  
  ZZZ_table<-merge(ZZZ_table,
                             ALL_dB_FINAL,
                             by=c("VAR"),
                             all=T)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_3\n")
    cat(str(ZZZ_table))
    cat("\n")
    cat(str(unique(ZZZ_table$VAR)))
    cat("\n")
  }
  
  
  
  ZZZ_table<-merge(ZZZ_table,
                             AF_and_CSQ_file_subset,
                             by=c("VAR","rs"),
                             all.x=T)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("ZZZ_table_4\n")
    cat(str(ZZZ_table))
    cat("\n")
    cat(str(unique(ZZZ_table$VAR)))
    cat("\n")
    cat(sprintf(as.character(colnames(ZZZ_table))))
    cat("\n")
  }
  
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/', sep='')
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  if(Condition_DEBUG == 1)
  {
    cat("path7:\t")
    cat(sprintf(as.character(path7)))
    cat("\n")
  }
  
  setwd(path7)
  
  write.table(ZZZ_table, sep="\t", quote=F, row.names = F,file="ZZZ_in_process.tsv")
  
  
  
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----


main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--CURATED_TABLE"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_INTERVAL_LM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_BP_LM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_LogLM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Proxy_R2_TABLE"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Proxy_CSQ_TABLE"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMMULATIVE_CLASSES"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_input"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genIE_active"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Prob_Threshold"), type="numeric", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--AF_and_CSQ_file"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
  ADD_INTERVAL_DE_Main_VARS(opt)
  ADD_INTERVAL_DTU_Main_VARS(opt)
  ADD_BP_DE_Main_VARS(opt)
  ADD_Proxys_CSQ(opt)
  ADD_Screenings(opt)
  ADD_phenotypes(opt)

}


###########################################################################

system.time( main() )
