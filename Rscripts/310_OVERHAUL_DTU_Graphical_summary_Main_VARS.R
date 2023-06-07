

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
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("gtools", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)


collate_results_function = function(option_list)
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
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS_UPDATED  = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS_UPDATED ))
  cat("\n")
  
  #### Read previous files ----
  
  setwd(out)
  
  Proxy_file<-as.data.frame(fread(file="Proxy_R2_TABLE_ALL.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("Proxy_file\n")
  cat(str(Proxy_file))
  cat("\n")
  cat(str(unique(Proxy_file$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file$Proxy_VAR)))
  cat("\n")
  
  
  ABSENT_WGS_RNA<-as.data.frame(fread(file="ABSENT_WGS_RNA_Seq.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("ABSENT_WGS_RNA\n")
  cat(str(ABSENT_WGS_RNA))
  cat("\n")
  
  
  #### UPDATE Selected_vars ----
  
  indx.UPDATE<-which(SELECTED_VARS_UPDATED %in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS_UPDATED [-indx.UPDATE]
    
  }else{
    
    SELECTED_VARS_UPDATED = SELECTED_VARS_UPDATED 
  }
  
  
  
  cat("SELECTED_VARS_UPDATED_\n")
  cat(str(SELECTED_VARS_UPDATED))
  cat("\n")
  
  
  
  Proxy_file_UPDATED<- Proxy_file[-which(Proxy_file$VAR%in%ABSENT_WGS_RNA$VAR |
                                           Proxy_file$Proxy_VAR%in%ABSENT_WGS_RNA$VAR ),]
  
  
  
  cat("Proxy_file_UPDATED\n")
  cat(str(Proxy_file_UPDATED))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$VAR)))
  cat("\n")
  cat(str(unique(Proxy_file_UPDATED$Proxy_VAR)))
  cat("\n")
  
  ######################### MASTER LOOP #############################
  
  Condition_DEBUG <- 0
  
  # SELECTED_VARS_UPDATED<-c("chr8_41589736_T_G","chr7_101499930_G_A","chr19_15653669_T_C")
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {

    
    
    #### Read transposed expression ----
    
    
    
    Transposed_Isoform_Expression<-readRDS(file=opt$Transposed_Isoform_Expression)
    colnames(Transposed_Isoform_Expression)[which(colnames(Transposed_Isoform_Expression) == "Sample_id")]<-"sample_id"
    
    #### Read Results_LogLM ----
    
    Results_LogLM<-as.data.frame(fread(file=opt$Results_LogLM, sep="\t", header=T) , stringsAsFactors=T)
    
    
    cat("Results_LogLM\n")
    cat(str(Results_LogLM))
    cat("\n")
    cat(str(unique(Results_LogLM$VAR)))
    cat("\n")
    
    ### size scale ----
    
    minuslospval_vector_DTU<-unique(c(Results_LogLM$CIS_gene_minuslogpvalue,Results_LogLM$Block_PCHiC_minuslogpvalue))
    
    cat("minuslospval_vector_DTU_0\n")
    cat(str(minuslospval_vector_DTU))
    cat("\n")
    
    minuslospval_vector_DTU<-minuslospval_vector_DTU[!is.na(minuslospval_vector_DTU)]
    
    cat("minuslospval_vector_DTU_1\n")
    cat(str(minuslospval_vector_DTU))
    cat("\n")
    
    minuslospval_vector_DTU_SIG<-minuslospval_vector_DTU[which(minuslospval_vector_DTU >= 1.3)]
    
    cat("minuslospval_vector_DTU_SIG_0\n")
    cat(str(minuslospval_vector_DTU_SIG))
    cat("\n")
    
    SUMMARY_minuslospval_vector_DTU_SIG<-summary(minuslospval_vector_DTU_SIG)
    
    cat("SUMMARY_minuslospval_vector_DTU_SIG\n")
    cat(sprintf(as.character(names(SUMMARY_minuslospval_vector_DTU_SIG))))
    cat("\n")
    cat(sprintf(as.character(SUMMARY_minuslospval_vector_DTU_SIG)))
    cat("\n")
    
    
    breaks.size<-sort(c(1,seq(0,8, by=2)))
    labels.size<-as.character(breaks.size)
    
    cat("labels.size\n")
    cat(str(labels.size))
    cat("\n")
    
    ##### LOOP TO READ ALL VARIABLES -----
    
   
    
    
    for(i in 1:length(SELECTED_VARS_UPDATED))
    {
      
      SELECTED_VARS_sel<-SELECTED_VARS_UPDATED[i]
      
      cat("------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(SELECTED_VARS_sel)))
      cat("\n")
      
      #### Read files for Residuals violin plots ----
      
      path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
      
      
    
      
      path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
      
      # cat("path7\n")
      # cat(sprintf(as.character(path7)))
      # cat("\n")
      
      
      if (file.exists(path7)){
        
        
        
        
      } else {
        dir.create(file.path(path7))
        
      }
      
      
      
      path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
      path8<-paste(path7,SELECTED_VARS_sel,'/', sep='')
      
      if (file.exists(path8)){
        
        
        
        
      } else {
        dir.create(file.path(path8))
        
      }
      
      setwd(path8)
      
      
      
      ### INTERVAL
      
      Results_LogLM_sel<-Results_LogLM[which(Results_LogLM$VAR%in%SELECTED_VARS_sel),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_LogLM_sel_\n")
        cat(str(Results_LogLM_sel))
        cat("\n")
        cat(str(unique(Results_LogLM_sel$VAR)))
        cat("\n")
      }
      
      
      ### ONLY  INTERVAL exists
      
      if(dim(Results_LogLM_sel)[1] >0)
      {
        ### Read Tappas_gff----
        
        
        Tappas_gff<-as.data.frame(fread(file=opt$Tappas_gff, sep="\t", header=F) , stringsAsFactors=F)
        
        colnames(Tappas_gff)<-c("transcript_id","annotation_source","type","START","END","DUMMY_1","strand","DUMMY2","Desc")
        
        
        # cat("Tappas_gff\n")
        # cat(str(Tappas_gff))
        # cat("\n")
        
        Tappas_gff_transcript<-Tappas_gff[which(Tappas_gff$type == "transcript"),]
        
        # cat("Tappas_gff_transcript_0\n")
        # cat(str(Tappas_gff_transcript))
        # cat("\n")
        
        Tappas_gff_transcript$BIOTYPE<-gsub("^.+primary_class=","",Tappas_gff_transcript$Desc)
        Tappas_gff_transcript$BIOTYPE<-gsub(";.+$","",Tappas_gff_transcript$BIOTYPE)
        
        indx.int<-c(which(colnames(Tappas_gff_transcript) == "transcript_id"),
                    which(colnames(Tappas_gff_transcript) == "BIOTYPE"))
        
        Tappas_gff_transcript<-unique(Tappas_gff_transcript[,indx.int])
        
        # cat("Tappas_gff_transcript_1\n")
        # cat(str(Tappas_gff_transcript))
        # cat("\n")
        
        Tappas_gff_CDS<-Tappas_gff[which(Tappas_gff$type == "CDS"),]
        
        # cat("Tappas_gff_CDS\n")
        # cat(str(Tappas_gff_CDS))
        # cat("\n")
        
        Tappas_gff_CDS$UniProt_ID<-gsub("^.+ID=","",Tappas_gff_CDS$Desc)
        Tappas_gff_CDS$UniProt_ID<-gsub(";.+$","",Tappas_gff_CDS$UniProt_ID)
        Tappas_gff_CDS$UniProt_ID<-gsub("ID=","ProtID=",Tappas_gff_CDS$UniProt_ID)
        
        
        indx.int<-c(which(colnames(Tappas_gff_CDS) == "transcript_id"),
                    which(colnames(Tappas_gff_CDS) == "UniProt_ID"))
        
        Tappas_gff_CDS<-unique(Tappas_gff_CDS[,indx.int])
        
        # cat("Tappas_gff_CDS_1\n")
        # cat(str(Tappas_gff_CDS))
        # cat("\n")
        
        Tappas_gff_DOMAIN<-Tappas_gff[which(Tappas_gff$type == "DOMAIN"),]
        
        
        # cat("Tappas_gff_DOMAIN\n")
        # cat(str(Tappas_gff_DOMAIN))
        # cat("\n")
        
        Tappas_gff_DOMAIN$PFAM_ID<-gsub("^.+ID=","",Tappas_gff_DOMAIN$Desc)
        Tappas_gff_DOMAIN$PFAM_ID<-gsub(";.+$","",Tappas_gff_DOMAIN$PFAM_ID)
        Tappas_gff_DOMAIN$PFAM_ID<-gsub("ID=","",Tappas_gff_DOMAIN$PFAM_ID)
        
        Tappas_gff_DOMAIN$PFAM_Desc<-gsub("^.+Desc=","",Tappas_gff_DOMAIN$Desc)
        Tappas_gff_DOMAIN$PFAM_Desc<-gsub(";.+$","",Tappas_gff_DOMAIN$PFAM_Desc)
        
        Tappas_gff_DOMAIN<-Tappas_gff_DOMAIN[order(Tappas_gff_DOMAIN$transcript_id,
                                                   Tappas_gff_DOMAIN$PFAM_ID),]
        
        Tappas_gff_DOMAIN$string_PFAM<-paste(Tappas_gff_DOMAIN$PFAM_ID,Tappas_gff_DOMAIN$PFAM_Desc,sep="__")
        
        indx.int<-c(which(colnames(Tappas_gff_DOMAIN) == "transcript_id"),
                    which(colnames(Tappas_gff_DOMAIN) == "string_PFAM"))
        
        
        
        Tappas_gff_DOMAIN<-unique(Tappas_gff_DOMAIN[,indx.int])
        
        
        # cat("Tappas_gff_DOMAIN_1\n")
        # cat(str(Tappas_gff_DOMAIN))
        # cat("\n")
        
        Tappas_gff_miRNA_Binding<-Tappas_gff[which(Tappas_gff$type == "miRNA_Binding"),]
        
        # cat("Tappas_gff_miRNA_Binding\n")
        # cat(str(Tappas_gff_miRNA_Binding))
        # cat("\n")
        Tappas_gff_miRNA_Binding$mRNA_feature_ID<-gsub("^.+ID=","",Tappas_gff_miRNA_Binding$Desc)
        Tappas_gff_miRNA_Binding$mRNA_feature_ID<-gsub(";.+$","",Tappas_gff_miRNA_Binding$mRNA_feature_ID)
        Tappas_gff_miRNA_Binding$mRNA_feature_ID<-gsub("ID=","",Tappas_gff_miRNA_Binding$mRNA_feature_ID)
        
        Tappas_gff_miRNA_Binding<-Tappas_gff_miRNA_Binding[order(Tappas_gff_miRNA_Binding$transcript_id,
                                                                 Tappas_gff_miRNA_Binding$mRNA_feature_ID),]
        
        indx.int<-c(which(colnames(Tappas_gff_miRNA_Binding) == "transcript_id"),
                    which(colnames(Tappas_gff_miRNA_Binding) == "mRNA_feature_ID"))
        
        
        
        Tappas_gff_miRNA_Binding<-unique(Tappas_gff_miRNA_Binding[,indx.int])
        
        
        # cat("Tappas_gff_miRNA_Binding_1\n")
        # cat(str(Tappas_gff_miRNA_Binding))
        # cat("\n")
        
        Tappas_gff_PTM<-Tappas_gff[which(Tappas_gff$type == "PTM"),]
        
        # cat("Tappas_gff_PTM\n")
        # cat(str(Tappas_gff_PTM))
        # cat("\n")
        
        Tappas_gff_PTM$PTM_residue<-gsub("^.+ID=","",Tappas_gff_PTM$Desc)
        Tappas_gff_PTM$PTM_residue<-gsub(";.+$","",Tappas_gff_PTM$PTM_residue)
        Tappas_gff_PTM$PTM_residue<-gsub("ID=","",Tappas_gff_PTM$PTM_residue)
        
        Tappas_gff_PTM$START<-as.numeric(Tappas_gff_PTM$START)
        
        Tappas_gff_PTM<-Tappas_gff_PTM[order(Tappas_gff_PTM$transcript_id,
                                             Tappas_gff_PTM$START,
                                             Tappas_gff_PTM$PTM_residue),]
        
        Tappas_gff_PTM$string_PTM<-paste(Tappas_gff_PTM$PTM_residue,Tappas_gff_PTM$START,Tappas_gff_PTM$END,sep="_")
        
        indx.int<-c(which(colnames(Tappas_gff_PTM) == "transcript_id"),
                    which(colnames(Tappas_gff_PTM) == "string_PTM"))
        
        
        
        Tappas_gff_PTM<-unique(Tappas_gff_PTM[,indx.int])
        
        # cat("Tappas_gff_PTM_1\n")
        # cat(str(Tappas_gff_PTM))
        # cat("\n")
        
        rm(Tappas_gff)
        
        
        setwd(path6)
        
        INTERVAL_covariates_and_PEER_factors_sel<-readRDS(file=paste("INTERVAL_covariates_and_PEER_factors_",SELECTED_VARS_sel,".rds", sep=''))
        
        # cat("INTERVAL_covariates_and_PEER_factors_sel\n")
        # cat(str(INTERVAL_covariates_and_PEER_factors_sel))
        # cat("\n")
        
        
        
        if (file.exists(paste("DTU_LogRatioLM_HET_RESIDUALS_",SELECTED_VARS_sel,".csv", sep=''))) {
          
          setwd(path6)
          
          residuals_df = as.data.frame(read.csv(file=paste("DTU_LogRatioLM_HET_RESIDUALS_",SELECTED_VARS_sel,".csv", sep=''), header=T), stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("residuals_df\n")
            cat(str(residuals_df))
            cat("\n")
          }
          
          residuals_df_sel<-residuals_df[which(residuals_df$transcript_id%in%Results_LogLM_sel$transcript_id),]
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("residuals_df_sel\n")
          #   cat(str(residuals_df_sel))
          #   cat("\n")
          # }
          
          
          

         
          ######################################################################
          
          ENSG_array<-unique(c(Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$CIS_gene_minuslogpvalue)],
                               Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$Block_PCHiC_minuslogpvalue)]))
          
          
          
          if(length(ENSG_array) >0)
          {
            cat("ENSG_array_\n")
            cat(str(ENSG_array))
            cat("\n")
            
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
              
              #"ajusted.minuslogpvalue_Genotypes","coefficient_Genotypes",
              
              df_aggregate<-data.frame()
              
              if(dim(CIS_gene_INTERVAL)[1] >0)
              {
                Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "transcript_id"),
                                                                     which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                     which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                     which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                                     which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"),
                                                                     which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta_Z_score"))]
                colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
                colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
                colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
                
                
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
                                                               which(colnames(CIS_gene_INTERVAL) == "CIS_Beta"),
                                                               which(colnames(CIS_gene_INTERVAL) == "CIS_Beta_Z_score"))]
                colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
                colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta")]<-"coefficient_Genotypes"
                colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
                
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
                                                                     which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"),
                                                                     which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta_Z_score"))]
                colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
                colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
                colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
                
                
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
                                                                                                                                                      which(colnames(Results_LogLM_sel_ENSG_sel) == "Block_PCHiC_Beta"),
                                                                                                                                                      which(colnames(Results_LogLM_sel_ENSG_sel) == "Block_PCHiC_Beta_Z_score"))]
              colnames(rescue_transcript_INTERVAL)[which(colnames(rescue_transcript_INTERVAL)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
              colnames(rescue_transcript_INTERVAL)[which(colnames(rescue_transcript_INTERVAL)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
              colnames(rescue_transcript_INTERVAL)[which(colnames(rescue_transcript_INTERVAL)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
              
              
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
              
              df_aggregate$VAR<-SELECTED_VARS_sel
              df_aggregate$Significance<-NA
              df_aggregate$RNASeq_source<-NA
              
              if(Condition_DEBUG == 1)
              {
                cat("df_aggregate_2\n")
                cat(str(df_aggregate))
                cat("\n")
                
              }
              
              #### Significance ----
              
              df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes >= 1.3)]<-"YES"
              df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes < 1.3)]<-"NO"
              
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
                cat("df_aggregate$RNASeq_source\n")
                cat(sprintf(as.character(names(summary(df_aggregate$RNASeq_source)))))
                cat("\n")
                cat(sprintf(as.character(summary(df_aggregate$RNASeq_source))))
                cat("\n")
              }
              
              ######## TPM, Ratio TPM and residuals tables -----
              
              Transposed_Isoform_Expression_sel<-Transposed_Isoform_Expression[,c(which(colnames(Transposed_Isoform_Expression) == "sample_id"),
                                                                                  which(colnames(Transposed_Isoform_Expression) %in% df_aggregate$transcript_id))]
              if(Condition_DEBUG == 1)
              {
                cat("Transposed_Isoform_Expression_sel_1\n")
                cat(str(Transposed_Isoform_Expression_sel))
                cat("\n")
              }
              
              ################# melt the Transposed_Isoform_Expression_sel & calculate Ratio ----
              
              Transposed_Isoform_Expression_sel.m<-melt(Transposed_Isoform_Expression_sel, id.vars=c("sample_id"),value.name = "log2TPM",variable.name = "transcript_id")
              Transposed_Isoform_Expression_sel.m$TPM<--1+2^(Transposed_Isoform_Expression_sel.m$log2TPM)
              
              
              if(Condition_DEBUG == 1)
              {
                cat("Transposed_Isoform_Expression_sel.m_0\n")
                cat(str(Transposed_Isoform_Expression_sel.m))
                cat("\n")
                
                cat("sample_id\n")
                cat(str(Transposed_Isoform_Expression_sel.m$sample_id))
                cat("\n")
                
                cat("transcript_id\n")
                cat(str(Transposed_Isoform_Expression_sel.m$transcript_id))
                cat("\n")
                
                # ########################################
                # quit(status = 1)
              }
              
              residuals_df_sel_ENSG_sel<-residuals_df_sel[which(residuals_df_sel$transcript_id%in%df_aggregate$transcript_id),]
              
              residuals_df_sel_ENSG_sel.m<-melt(residuals_df_sel_ENSG_sel, id.vars=c("transcript_id"),value.name = "residuals_Ratio",variable.name = "sample_id")
              
              
              if(Condition_DEBUG == 1)
              {
                # cat("residuals_df_sel_ENSG_sel\n")
                # cat(str(residuals_df_sel_ENSG_sel))
                # cat("\n")
                
                cat("residuals_df_sel_ENSG_sel.m\n")
                cat(str(residuals_df_sel_ENSG_sel.m))
                cat("\n")
              }
              
              
              
              Transposed_Isoform_Expression_sel.m<-merge(Transposed_Isoform_Expression_sel.m,
                                                         residuals_df_sel_ENSG_sel.m,
                                                         by=c("sample_id","transcript_id"),
                                                         all=T)
              
              
              if(Condition_DEBUG == 1)
              {
                cat("Transposed_Isoform_Expression_sel.m_1\n")
                cat(str(Transposed_Isoform_Expression_sel.m))
                cat("\n")
                
                # #####################################################
                # quit(status = 1)
              }
              
              A<-round(summary(Transposed_Isoform_Expression_sel.m$residuals_Ratio[!is.na(Transposed_Isoform_Expression_sel.m$residuals_Ratio)]),2)
              # A2<-round(summary(Transposed_Isoform_Expression_sel.m$TPM[!is.na(Transposed_Isoform_Expression_sel.m$TPM)]),2)
              
              Transposed_Isoform_Expression_sel.m$residuals_Ratio_no_negative<-Transposed_Isoform_Expression_sel.m$residuals_Ratio+abs(A[1])
              
              A3<-round(summary(Transposed_Isoform_Expression_sel.m$residuals_Ratio_no_negative[!is.na(Transposed_Isoform_Expression_sel.m$residuals_Ratio_no_negative)]),0)
              
              
              
              if(Condition_DEBUG == 1)
              {
                cat("Transposed_Isoform_Expression_sel.m_2\n")
                cat(str(Transposed_Isoform_Expression_sel.m))
                cat("\n")
                
                cat("Summary_residuals_Ratio:\t")
                cat(sprintf(as.character(names(A))))
                cat("\n")
                cat(sprintf(as.character(A)))
                cat("\n")
                
                # cat("Summary_TPM:\t")
                # cat(sprintf(as.character(names(A2))))
                # cat("\n")
                # cat(sprintf(as.character(A2)))
                # cat("\n")
                
                
                cat("Summary_residuals_Ratio_no_negative:\t")
                cat(sprintf(as.character(names(A3))))
                cat("\n")
                cat(sprintf(as.character(A3)))
                cat("\n")
                
                
              }
              
              
              
              ### calculate summatory TPM per sample
              
              Transposed_Isoform_Expression_sel.m.dt<-data.table(Transposed_Isoform_Expression_sel.m,
                                                                 key=c("sample_id"))
              
              if(Condition_DEBUG == 1)
              {
                cat("Transposed_Isoform_Expression_sel.m.dt_0\n")
                cat(str(Transposed_Isoform_Expression_sel.m.dt))
                cat("\n")
              }
              
              
              Summary_table_GENE_EXP<-as.data.frame(Transposed_Isoform_Expression_sel.m.dt[, .(sum_GENE_EXP=sum(TPM)),
                                                                                           by=key(Transposed_Isoform_Expression_sel.m.dt)],stringsAsFactors=F)
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Summary_table_GENE_EXP_0\n")
                cat(str(Summary_table_GENE_EXP))
                cat("\n")
                
                # quit(status = 1)
              }
              
              
              Transposed_Isoform_Expression_sel.m<-merge(Transposed_Isoform_Expression_sel.m,
                                                         Summary_table_GENE_EXP,
                                                         by="sample_id")
              
              
              if(Condition_DEBUG == 1)
              {
                cat("Transposed_Isoform_Expression_sel.m_1\n")
                cat(str(Transposed_Isoform_Expression_sel.m))
                cat("\n")
                
              }
              
              Transposed_Isoform_Expression_sel.m.dt<-data.table(Transposed_Isoform_Expression_sel.m,
                                                                 key=c("sample_id","transcript_id"))
              
              Ratio_df<-as.data.frame(Transposed_Isoform_Expression_sel.m.dt[, .(log2TPM=log2TPM,
                                                                                 TPM=TPM,
                                                                                 sum_GENE_EXP=sum_GENE_EXP,
                                                                                 residuals_Ratio_no_negative=residuals_Ratio_no_negative,
                                                                                 Ratio=TPM/sum_GENE_EXP),by=key(Transposed_Isoform_Expression_sel.m.dt)],stringsAsFactors=F)
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_0\n")
                cat(str(Ratio_df))
                cat("\n")
                # 
                # cat("distrib_ratios\n")
                # cat(sprintf(as.character(names(summary(Ratio_df$Ratio)))))
                # cat("\n")
                # cat(sprintf(as.character(summary(Ratio_df$Ratio))))
                # cat("\n")
                # 
                # 
                # quit(status = 1)
              }
              
              
              
              
              
              #### Merge with covariates matrix & keep HET ----
              
              
              
              Ratio_df<-merge(Ratio_df,
                              INTERVAL_covariates_and_PEER_factors_sel,
                              by=c("sample_id"),
                              all=T)
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_3\n")
                cat(str(Ratio_df))
                cat("\n")
                
                # quit(status = 1)
                
              }
              
              Ratio_df_HET<-droplevels(Ratio_df[which(Ratio_df$Genotype != "HOM"),])
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_HET_0\n")
                cat(str(Ratio_df_HET))
                cat("\n")
                
                
                # quit(status = 1)
                
              }
              
              
              ##### Violin visual log2TPM ----
              
              A<-round(summary(Ratio_df_HET$log2TPM[!is.na(Ratio_df_HET$log2TPM)]),2)
              
              
              # cat("summary_GENE_EXP\n")
              # cat(sprintf(as.character(names(A))))
              # cat("\n")
              # cat(sprintf(as.character(A)))
              # cat("\n")
              
              step<-abs(A[6]-A[1])/10
              
              if(step == 0)
              {
                
                step<-1
              }
              
              breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
              labels.Rank<-as.character(round(breaks.Rank,3))
              
              levels_transcript_id<-unique(Ratio_df_HET$transcript_id)
              
              if(Condition_DEBUG == 1)
              {
                cat("labels.Rank:\t")
                cat(sprintf(as.character(labels.Rank)))
                cat("\n")
                
                cat("levels_transcript_id:\t")
                cat(sprintf(as.character(levels_transcript_id)))
                cat("\n")
                
                # quit(status = 1)
              }
              
              
              
              
              graph_log2TPM<-ggplot(data=Ratio_df_HET,aes(x=Genotype, y=log2TPM, fill=Genotype)) +
                geom_violin()+
                stat_summary(fun = median, fun.min = median, fun.max = median,
                             geom = "crossbar", width = 0.5)+
                theme_bw()+
                theme(axis.title.y=element_text(size=24, family="sans"),
                      axis.title.x=element_text(size=24, family="sans"),
                      axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                      axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                      legend.title=element_text(size=16,color="black", family="sans"),
                      legend.text=element_text(size=12,color="black", family="sans"))+
                scale_x_discrete(name=NULL, drop=F)+
                scale_y_continuous(name="Transcript EXP normalised (log2TPM)",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
                scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
                ggtitle(paste(unique(Ratio_df_HET$ensembl_gene_id),unique(Ratio_df_HET$HGNC),sep=' '))+
                ggeasy::easy_center_title()
              
              # cat("facet_wrap_1:\t")
              # cat("\n")
              
              
              
              graph_log2TPM <- graph_log2TPM +
                facet_grid(cols = vars(Ratio_df_HET$transcript_id), scales = "free") +
                theme_cowplot(font_size = 16)+
                theme( strip.background = element_blank(),
                       strip.placement = "inside",
                       strip.text = element_text(size=15),
                       panel.spacing = unit(0.2, "lines"), 
                       panel.background=element_rect(fill="white"),
                       panel.border=element_rect(colour="black",size=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
                theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=12))
              
              ##### Violin visual Ratio ----
              
              # A<-round(summary(Ratio_df_HET$logRatio[!is.na(Ratio_df_HET$logRatio)]),5)
              # 
              # 
              # 
              # step<-abs(A[6]-A[1])/10
              # 
              # if(step == 0)
              # {
              #   
              #   step<-0.1
              # }
              
              breaks.Rank<-unique(seq(from= 0, to=1,by=0.1))
              labels.Rank<-as.character(breaks.Rank)
              
              levels_transcript_id<-unique(Ratio_df_HET$transcript_id)
              
              if(Condition_DEBUG == 1)
              {
                cat("labels.Rank:\t")
                cat(sprintf(as.character(labels.Rank)))
                cat("\n")
                
                cat("levels_transcript_id:\t")
                cat(sprintf(as.character(levels_transcript_id)))
                cat("\n")
                
                # quit(status = 1)
              }
              
              
              
              
              graph_Ratio<-ggplot(data=Ratio_df_HET,aes(x=Genotype, y=Ratio, fill=Genotype)) +
                geom_violin()+
                stat_summary(fun = median, fun.min = median, fun.max = median,
                             geom = "crossbar", width = 0.5)+
                theme_bw()+
                theme(axis.title.y=element_text(size=24, family="sans"),
                      axis.title.x=element_text(size=24, family="sans"),
                      axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                      axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                      legend.title=element_text(size=16,color="black", family="sans"),
                      legend.text=element_text(size=12,color="black", family="sans"))+
                scale_x_discrete(name=NULL, drop=F)+
                scale_y_continuous(name="Ratio Transcript TPM/geneTPM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
                scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
                ggtitle(paste(unique(Ratio_df_HET$ensembl_gene_id),unique(Ratio_df_HET$HGNC),sep=' '))+
                ggeasy::easy_center_title()
              
              # cat("facet_wrap_1:\t")
              # cat("\n")
              
              
              
              graph_Ratio <- graph_Ratio +
                facet_grid(cols = vars(Ratio_df_HET$transcript_id), scales = "free") +
                theme_cowplot(font_size = 16)+
                theme( strip.background = element_blank(),
                       strip.placement = "inside",
                       strip.text = element_text(size=15),
                       panel.spacing = unit(0.2, "lines"), 
                       panel.background=element_rect(fill="white"),
                       panel.border=element_rect(colour="black",size=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
                theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=12))
              
              
              ##### Violin visual Ratio_adjusted reduced model ----
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_HET_TPM_reduced_model\n")
                cat(str(Ratio_df_HET))
                cat("\n")
                
                
                # quit(status = 1)
                
              }
              
              A<-round(summary(Ratio_df_HET$residuals_Ratio_no_negative[!is.na(Ratio_df_HET$residuals_Ratio_no_negative)]),2)
              
              
              
              step<-abs(A[6]-A[1])/10
              
              if(step == 0)
              {
                
                step<-0.1
              }
              
              breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
              
              A<-summary(Ratio_df_HET$residuals_Ratio_no_negative)
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Summary\n")
                cat(sprintf(as.character(names(A))))
                cat("\n")
                cat(sprintf(as.character(A)))
                cat("\n")
                
                cat("breaks.Rank1\n")
                cat(sprintf(as.character(names(breaks.Rank))))
                cat("\n")
                cat(sprintf(as.character(breaks.Rank)))
                cat("\n")
                
                # quit(status = 1)
                
              }
              
              
              # breaks.Rank<-unique(seq(from= 0, to=1,by=0.1))
              labels.Rank<-as.character(breaks.Rank)
              
              levels_transcript_id<-unique(Ratio_df_HET$transcript_id)
              
              if(Condition_DEBUG == 1)
              {
                cat("labels.Rank2:\t")
                cat(sprintf(as.character(labels.Rank)))
                cat("\n")
                
                cat("levels_transcript_id:\t")
                cat(sprintf(as.character(levels_transcript_id)))
                cat("\n")
                
                # quit(status = 1)
              }
              
              
              
              
              graph_adjusted_Ratio<-ggplot(data=Ratio_df_HET,aes(x=Genotype, y=residuals_Ratio_no_negative, fill=Genotype)) +
                geom_violin()+
                stat_summary(fun = median, fun.min = median, fun.max = median,
                             geom = "crossbar", width = 0.5)+
                theme_bw()+
                theme(axis.title.y=element_text(size=24, family="sans"),
                      axis.title.x=element_text(size=24, family="sans"),
                      axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                      axis.text.x=element_text(angle=0,size=18, color="black", family="sans"),
                      legend.title=element_text(size=16,color="black", family="sans"),
                      legend.text=element_text(size=12,color="black", family="sans"))+
                scale_x_discrete(name=NULL, drop=F)+
                scale_y_continuous(name="Residuals Ratio Transcript TPM/geneTPM adjusted for covariates",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
                scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
                ggtitle(paste(unique(Ratio_df_HET$ensembl_gene_id),unique(Ratio_df_HET$HGNC),sep=' '))+
                ggeasy::easy_center_title()
              
              # cat("facet_wrap_1:\t")
              # cat("\n")
              
              
              
              graph_adjusted_Ratio <- graph_adjusted_Ratio +
                facet_grid(cols = vars(Ratio_df_HET$transcript_id), scales = "free") +
                theme_cowplot(font_size = 14)+
                theme( strip.background = element_blank(),
                       strip.placement = "inside",
                       strip.text = element_text(size=14),
                       panel.spacing = unit(0.2, "lines"), 
                       panel.background=element_rect(fill="white"),
                       panel.border=element_rect(colour="black",size=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
                theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=12))
              
              ####### Dot plot genes ----
              
              ### update break.size
              
              local_minuslogpval_max<-max(df_aggregate$ajusted.minuslogpvalue_Genotypes[!is.na(df_aggregate$ajusted.minuslogpvalue_Genotypes)])
              
              if(Condition_DEBUG == 1)
              {
                cat("local_minuslogpval_max\n")
                cat(str(local_minuslogpval_max))
                cat("\n")
              }
              
              breaks.size_bigger_than_local_max<-breaks.size[which(breaks.size > local_minuslogpval_max)]
              
              if(Condition_DEBUG == 1)
              {
                cat("breaks.size_bigger_than_local_max\n")
                cat(str(breaks.size_bigger_than_local_max))
                cat("\n")
              }
              
              if(length(breaks.size_bigger_than_local_max) >0)
              {
                breaks.size_updated<-breaks.size
                
              }else{
                
                breaks.size_updated<-sort(unique(c(breaks.size[-length(breaks.size)],local_minuslogpval_max+0.1)))
                
                cat("breaks.size_updated\n")
                cat(sprintf(as.character(breaks.size_updated)))
                cat("\n")
                
                # quit(status = 1)
              }
              
              labels.size_updated<-as.character(breaks.size_updated)
              
              
              ### Use Beta from the full model of logRatio ---
              
              ### Labels beta
              
              indx.finite<-is.finite(df_aggregate$coefficient_Genotypes_Z_score)
              
              check.finite<-sum(indx.finite)
              
              if(check.finite >0)
              {
                A<-summary(df_aggregate$coefficient_Genotypes_Z_score[indx.finite])
                
                min_value<-A[1]
                max_value<-A[6]
                
                if(Condition_DEBUG == 1)
                {
                  cat("Values\n")
                  cat(sprintf(as.character((min_value))))
                  cat("\n")
                  cat(sprintf(as.character((max_value))))
                  cat("\n")
                  
                  cat("Step\n")
                }
                
                
                step<-(max_value-min_value)/2
                
                
                if(Condition_DEBUG == 1)
                {
                  cat(sprintf(as.character((step))))
                  cat("\n")
                }
                
                candidate_vector<-seq(min_value,max_value+step, by=step)
                
                breaks.Beta_INTERVAL<-sort(unique(round(c(0,candidate_vector,-0.025,0.025),3)))
                labels.Beta_INTERVAL<-as.character(breaks.Beta_INTERVAL)
                
              }else{
                
                breaks.Beta_INTERVAL<-sort(unique(round(c(0,-0.025,0.025),3)))
                labels.Beta_INTERVAL<-as.character(breaks.Beta_INTERVAL)
                
              }# length(indx.finite) >0
              
              if(Condition_DEBUG == 1)
              {
                cat("breaks.Beta_INTERVAL\n")
                cat(str(breaks.Beta_INTERVAL))
                cat("\n")
                
                cat("labels.Beta_INTERVAL\n")
                cat(str(labels.Beta_INTERVAL))
                cat("\n")
                
                # quit(status = 1)
              }
              
              ### Tappas ----
              
              df_aggregate<-merge(df_aggregate,
                                  Tappas_gff_transcript,
                                  by="transcript_id",
                                  all.x=T)
              
              df_aggregate$transcript_id<-factor(df_aggregate$transcript_id,
                                                 levels=rev(levels_transcript_id),
                                                 ordered=T)
              
              df_aggregate$Significance<-"NO"
              df_aggregate$Significance[which(df_aggregate$ajusted.minuslogpvalue_Genotypes >= 1.3)]<-"YES"
              
              df_aggregate$Significance<-factor(df_aggregate$Significance,
                                                levels=c("NO","YES"),
                                                ordered=T)
              
              if(Condition_DEBUG == 1)
              {
                cat("df_aggregate_for_dot_plot_1:\t")
                cat(str(df_aggregate))
                cat("\n")
                
              }
              
              df_aggregate<-merge(df_aggregate,
                                  Tappas_gff_CDS,
                                  by="transcript_id",
                                  all.x=T)
              
              Tappas_gff_miRNA_Binding_sel<-Tappas_gff_miRNA_Binding[which(Tappas_gff_miRNA_Binding$transcript_id%in%df_aggregate$transcript_id),]
              
              if(Condition_DEBUG == 1)
              {
                # cat("Tappas_gff_miRNA_Binding_sel_0:\t")
                # cat(str(Tappas_gff_miRNA_Binding_sel))
                # cat("\n")
                
              }
              
              Tappas_gff_miRNA_Binding_sel.dt<-data.table(Tappas_gff_miRNA_Binding_sel, key="transcript_id")
              
              
              Tappas_gff_miRNA_Binding_sel_DEF<-unique(as.data.frame(Tappas_gff_miRNA_Binding_sel.dt[,.(string_mRNA_feature_ID=paste(mRNA_feature_ID, collapse=";")),
                                                                                                     by=key(Tappas_gff_miRNA_Binding_sel.dt)], stringsAsFactors=F))
              
              if(Condition_DEBUG == 1)
              {
                cat("Tappas_gff_miRNA_Binding_sel_DEF_0:\t")
                cat(str(Tappas_gff_miRNA_Binding_sel_DEF))
                cat("\n")
                
              }
              
              Tappas_gff_DOMAIN_sel<-Tappas_gff_DOMAIN[which(Tappas_gff_DOMAIN$transcript_id%in%df_aggregate$transcript_id),]
              
              if(Condition_DEBUG == 1)
              {
                # cat("Tappas_gff_DOMAIN_sel_0:\t")
                # cat(str(Tappas_gff_DOMAIN_sel))
                # cat("\n")
                
              }
              
              Tappas_gff_DOMAIN_sel.dt<-data.table(Tappas_gff_DOMAIN_sel, key="transcript_id")
              
              
              Tappas_gff_DOMAIN_sel_DEF<-unique(as.data.frame(Tappas_gff_DOMAIN_sel.dt[,.(string_PFAM=paste(string_PFAM, collapse=";")),
                                                                                       by=key(Tappas_gff_DOMAIN_sel.dt)], stringsAsFactors=F))
              
              if(Condition_DEBUG == 1)
              {
                cat("Tappas_gff_DOMAIN_sel_DEF_0:\t")
                cat(str(Tappas_gff_DOMAIN_sel_DEF))
                cat("\n")
                
              }
              
              Tappas_gff_PTM_sel<-Tappas_gff_PTM[which(Tappas_gff_PTM$transcript_id%in%df_aggregate$transcript_id),]
              
              if(Condition_DEBUG == 1)
              {
                # cat("Tappas_gff_PTM_sel_0:\t")
                # cat(str(Tappas_gff_PTM_sel))
                # cat("\n")
                
                # quit(status = 1)
              }
              
              Tappas_gff_PTM_sel.dt<-data.table(Tappas_gff_PTM_sel, key="transcript_id")
              
              
              Tappas_gff_PTM_sel_DEF<-unique(as.data.frame(Tappas_gff_PTM_sel.dt[,.(string_PTM=paste(string_PTM, collapse=";")),
                                                                                 by=key(Tappas_gff_PTM_sel.dt)], stringsAsFactors=F))
              
              if(Condition_DEBUG == 1)
              {
                cat("Tappas_gff_PTM_sel_DEF_0:\t")
                cat(str(Tappas_gff_PTM_sel_DEF))
                cat("\n")
                
              }
              
              Tappas_gff_DEF<-merge(Tappas_gff_DOMAIN_sel_DEF,
                                    Tappas_gff_PTM_sel_DEF,
                                    by=c("transcript_id"),
                                    all=T)
              
              if(Condition_DEBUG == 1)
              {
                cat("Tappas_gff_DEF_0:\t")
                cat(str(Tappas_gff_DEF))
                cat("\n")
                
                # quit(status = 1)
              }
              
              Tappas_gff_DEF<-merge(Tappas_gff_DEF,
                                    Tappas_gff_miRNA_Binding_sel_DEF,
                                    by=c("transcript_id"),
                                    all=T)
              
              if(Condition_DEBUG == 1)
              {
                cat("Tappas_gff_DEF_1:\t")
                cat(str(Tappas_gff_DEF))
                cat("\n")
                
                # quit(status = 1)
              }
              
              
              ### Graph ----
              
              df_aggregate<-df_aggregate[order(df_aggregate$transcript_id),]
              
              if(Condition_DEBUG == 1)
              {
                cat("df_aggregate_for_dot_plot_2:\t")
                cat(str(df_aggregate))
                cat("\n")
                
              }

              # geom_point(data=subset(df_aggregate, df_aggregate$Significance == "YES"),aes(color=Significance,fill=coefficient_Genotypes_Z_score,size=ajusted.minuslogpvalue_Genotypes), stroke=2, shape=21)+
                
              dotplot_INTERVAL_DTU<-df_aggregate %>%
                mutate(myaxis = paste0(transcript_id, "\n", BIOTYPE,"\n",UniProt_ID)) %>%
                mutate(myaxis=fct_reorder(myaxis,as.numeric(df_aggregate$transcript_id))) %>%
                ggplot(aes(y=myaxis,
                           x=RNASeq_source)) +
                geom_point(aes(color=Significance,fill=coefficient_Genotypes_Z_score,size=ajusted.minuslogpvalue_Genotypes), stroke=1, shape=21)+
                scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
                scale_size(range = c(0,20), name='-log10pval',
                           breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)])) +
                scale_fill_gradient2(
                  low = "blue", 
                  mid = "white", 
                  high = "red", 
                  midpoint = 0,
                  breaks=breaks.Beta_INTERVAL,labels=labels.Beta_INTERVAL,
                  limits=c(breaks.Beta_INTERVAL[1]-0.01,breaks.Beta_INTERVAL[length(breaks.Beta_INTERVAL)]+0.01),name=paste('Effect size',
                                                                                                                            'Z score',
                                                                                                                            sep="\n"),na.value = "gray")+
                scale_y_discrete(name=NULL, drop=F)+
                scale_x_discrete(name=NULL, drop=F)+
                theme_classic()+
                ggeasy::easy_center_title()
              
              
              dotplot_INTERVAL_DTU<-dotplot_INTERVAL_DTU+
                facet_grid(cols = vars(df_aggregate$RNASeq_source), scales='free_x', space='free_x') +
                theme_cowplot(font_size = 14)+
                theme( strip.background = element_blank(),
                       strip.placement = "inside",
                       strip.text = element_text(size=14),
                       panel.spacing = unit(0.2, "lines"),
                       panel.background=element_rect(fill="white"),
                       panel.border=element_rect(colour="black",size=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
                theme(axis.title.y=element_blank(),
                      axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())+
                theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
                      legend.key.height = unit(1.5, 'cm'), #change legend key height
                      legend.key.width = unit(1, 'cm'), #change legend key width
                      legend.title = element_text(size=14), #change legend title font size
                      legend.text = element_text(size=14))+ #change legend text font size
                scale_x_discrete(name=NULL, drop=T)
              
              
              
              ### PRINT GRAPH ----
              path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
              path8<-paste(path7,SELECTED_VARS_sel,'/', sep='')
              path9<-paste(path8,HGNC_sel,'/', sep='')
              
              if (file.exists(path9)){
                
                
                
                
              } else {
                dir.create(file.path(path9))
                
              }
              
              setwd(path9)
              
              graph_DEF<-plot_grid(graph_log2TPM,graph_Ratio,graph_adjusted_Ratio,dotplot_INTERVAL_DTU,
                                   nrow = 2,
                                   ncol = 2,
                                   rel_widths=c(1,1),
                                   rel_heights=c(1,1))
              
             
              svgname<-paste("DTU_EVALUATION_",HGNC_sel,".svg",sep='')
              makesvg = TRUE
              
              if (makesvg == TRUE)
              {
                ggsave(svgname, plot= graph_DEF,
                       device="svg",
                       height=10, width=12)
              }
              
              graph_adjusted_Ratio <- graph_adjusted_Ratio +
                theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=14))+
                theme(axis.text.x=element_blank())
              
              
              graph_DEF<-plot_grid(graph_adjusted_Ratio,dotplot_INTERVAL_DTU,
                                   nrow = 1,
                                   ncol = 2,
                                   rel_widths=c(1,0.75))
              
              
              svgname<-paste("DTU_readout_",HGNC_sel,".svg",sep='')
              makesvg = TRUE
              
              if (makesvg == TRUE)
              {
                ggsave(svgname, plot= graph_DEF,
                       device="svg",
                       height=10, width=12)
              }
              
              
              
              
              if(Condition_DEBUG == 1)
              {
                cat("PRINTING#1\n")
                
              }
              # quit(status = 1)
              
              
           
              write.table(Tappas_gff_DEF, file=paste("Tappas_annotation_per_transcript_",HGNC_sel,".tsv",sep=''), sep="\t",quote=F, row.names = F)
              
              if(Condition_DEBUG == 1)
              {
                
                cat("SAVING:\n")
                
                # quit(status = 1)
              }
              
              # quit(status = 1)
            }# z 1:length(ENSG_array)
          }else{
            
            # no gene in CIS or Block
            
          }# gene in CIS or Block
          
          
        }#file.exists(paste("DTU_LogRatioLM_HET_RESIDUALS_",SELECTED_VARS_sel,".csv", sep=''))
      }# INTERVAL exist
      
    }# i in 1:length(SELECTED_VARS)
    
  }#length(SELECTED_VARS_UPDATED) >0
 
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
    make_option(c("--SELECTED_VARS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_LogLM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Transposed_Isoform_Expression"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Tappas_gff"), type="character", default=NULL,
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
  
 
  collate_results_function(opt)
 
 
  
  
  
}


###########################################################################

system.time( main() )
