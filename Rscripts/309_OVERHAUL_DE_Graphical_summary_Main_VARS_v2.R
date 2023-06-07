

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

graphical_function = function(option_list)
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
  
  #### READ and transform path_BP ----
  
  path_BP = opt$path_BP
  
  cat("path_BP\n")
  cat(sprintf(as.character(path_BP)))
  cat("\n")
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS_UPDATED  = unlist(strsplit(opt$SELECTED_VARS_UPDATED , split=","))
  
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
  
 
  
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    #### Read Thomas matrix of gene expression ----
    
    INTERVAL_GENE_EXP<-as.data.frame(fread(file=opt$INTERVAL_GENE_EXP, sep=",", header=T) , stringsAsFactors=F)
    
    
    
    # cat("INTERVAL_GENE_EXP\n")
    # cat(str(INTERVAL_GENE_EXP))
    # cat("\n")
    # 
    Samples<-unique(INTERVAL_GENE_EXP$sample_id)
    
    # cat("Samples\n")
    # cat(str(Samples))
    # cat("\n")
    
    ENSG_array<-unique(colnames(INTERVAL_GENE_EXP)[-which(colnames(INTERVAL_GENE_EXP) == "sample_id")])
    
    cat("ENSG_array\n")
    cat(str(ENSG_array))
    cat("\n")
    
    
    
   
    
    
    # quit(status=1)
    
    
    ### Read Results_INTERVAL_LM----
    
    
    
    Results_INTERVAL_LM<-as.data.frame(fread(file=opt$Results_INTERVAL_LM, sep="\t", header=T) , stringsAsFactors=T)
    
    
    cat("Results_INTERVAL_LM\n")
    cat(str(Results_INTERVAL_LM))
    cat("\n")
    cat(str(unique(Results_INTERVAL_LM$VAR)))
    cat("\n")
    
    ### Read Results_BP_LM----
    
    
    
    Results_BP_LM<-as.data.frame(fread(file=opt$Results_BP_LM, sep="\t", header=T) , stringsAsFactors=T)
    
    
    cat("Results_BP_LM\n")
    cat(str(Results_BP_LM))
    cat("\n")
    cat(str(unique(Results_BP_LM$VAR)))
    cat("\n")
    
    
   
    
    
    ### size scale ----
    
    minuslospval_vector_INTERVAL<-unique(c(Results_INTERVAL_LM$CIS_gene_minuslogpvalue,Results_INTERVAL_LM$Block_PCHiC_minuslogpvalue))
    
    cat("minuslospval_vector_INTERVAL_0\n")
    cat(str(minuslospval_vector_INTERVAL))
    cat("\n")
    
    minuslospval_vector_INTERVAL<-minuslospval_vector_INTERVAL[!is.na(minuslospval_vector_INTERVAL)]
    
    cat("minuslospval_vector_INTERVAL_1\n")
    cat(str(minuslospval_vector_INTERVAL))
    cat("\n")
    
    minuslospval_vector_INTERVAL_SIG<-minuslospval_vector_INTERVAL[which(minuslospval_vector_INTERVAL >= 1.3)]
    
    cat("minuslospval_vector_INTERVAL_SIG_0\n")
    cat(str(minuslospval_vector_INTERVAL_SIG))
    cat("\n")
    
    SUMMARY_minuslospval_vector_INTERVAL_SIG<-summary(minuslospval_vector_INTERVAL_SIG)
    
    cat("SUMMARY_minuslospval_vector_INTERVAL_SIG\n")
    cat(sprintf(as.character(names(SUMMARY_minuslospval_vector_INTERVAL_SIG))))
    cat("\n")
    cat(sprintf(as.character(SUMMARY_minuslospval_vector_INTERVAL_SIG)))
    cat("\n")
    
    minuslospval_vector_BP<-unique(c(Results_BP_LM$CIS_gene_minuslogpvalue,Results_BP_LM$Block_PCHiC_minuslogpvalue))
    
    cat("minuslospval_vector_BP_0\n")
    cat(str(minuslospval_vector_BP))
    cat("\n")
    
    minuslospval_vector_BP<-minuslospval_vector_BP[!is.na(minuslospval_vector_BP)]
    
    cat("minuslospval_vector_BP_1\n")
    cat(str(minuslospval_vector_BP))
    cat("\n")
    
    
    minuslospval_vector_BP_SIG<-minuslospval_vector_BP[which(minuslospval_vector_BP >= 1.3)]
    
    cat("minuslospval_vector_BP_SIG_0\n")
    cat(str(minuslospval_vector_BP_SIG))
    cat("\n")
    
    SUMMARY_minuslospval_vector_BP_SIG<-summary(minuslospval_vector_BP_SIG)
    
    cat("SUMMARY_minuslospval_vector_BP_SIG\n")
    cat(sprintf(as.character(names(SUMMARY_minuslospval_vector_BP_SIG))))
    cat("\n")
    cat(sprintf(as.character(SUMMARY_minuslospval_vector_BP_SIG)))
    cat("\n")
    
    breaks.size<-sort(c(1,seq(0,8, by=2)))
    labels.size<-as.character(breaks.size)
    
    cat("labels.size\n")
    cat(str(labels.size))
    cat("\n")
    
    ##### LOOP TO READ ALL VARIABLES -----
    
    Condition_DEBUG <- 0
    
    # SELECTED_VARS_UPDATED<-c("chr17_38764524_T_A","chr3_71355240_G_C","chr12_111844956_C_T")
    
    # SELECTED_VARS_UPDATED<-"chr19_15653669_T_C"
    
    
    list<-list()
    List_INTERVAL<-list()
    
    
    
    for(i in 1:length(SELECTED_VARS_UPDATED))
    {
      
      SELECTED_VARS_sel<-SELECTED_VARS_UPDATED[i]
      
      cat("------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(SELECTED_VARS_sel)))
      cat("\n")
      
      #### Read files for Residuals violin plots ----
      
      path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
      
      
      setwd(path6)
      
      # cat("Path_and_files\n")
      # cat(sprintf(as.character(path6)))
      # cat("\n")
      # cat(sprintf(as.character(paste("DE_LM_HET_RESIDUALS_",SELECTED_VARS_sel,".csv", sep=''))))
      # cat("\n")
      
      INTERVAL_covariates_and_PEER_factors_sel<-readRDS(file=paste("INTERVAL_covariates_and_PEER_factors_",SELECTED_VARS_sel,".rds", sep=''))
      
      # cat("INTERVAL_covariates_and_PEER_factors_sel\n")
      # cat(str(INTERVAL_covariates_and_PEER_factors_sel))
      # cat("\n")
     
      
      residuals_INTERVAL_df = read.csv(file=paste("DE_LM_HET_RESIDUALS_",SELECTED_VARS_sel,".csv", sep=''),sep=",", header=T)
      
      # cat("residuals_INTERVAL_df_0\n")
      # cat(str(residuals_INTERVAL_df))
      # cat("\n")
      
      # quit(status = 1)
      
      
      path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
      
      # cat("path7\n")
      # cat(sprintf(as.character(path7)))
      # cat("\n")
      
      
      if (file.exists(path7)){
        
        
        
        
      } else {
        dir.create(file.path(path7))
        
      }
      
      
      
      
      path8<-paste(path7,SELECTED_VARS_sel,'/', sep='')
      
      if (file.exists(path8)){
        
        
        
        
      } else {
        dir.create(file.path(path8))
        
      }
      
      setwd(path8)
      
      
      ### INTERVAL
      
      Results_INTERVAL_LM_sel<-Results_INTERVAL_LM[which(Results_INTERVAL_LM$VAR%in%SELECTED_VARS_sel),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_INTERVAL_LM_sel_\n")
        cat(str(Results_INTERVAL_LM_sel))
        cat("\n")
        cat(str(unique(Results_INTERVAL_LM_sel$VAR)))
        cat("\n")
      }
      
      ### BP
      
      Results_BP_LM_sel<-Results_BP_LM[which(Results_BP_LM$VAR%in%SELECTED_VARS_sel),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_BP_LM_sel_\n")
        cat(str(Results_BP_LM_sel))
        cat("\n")
        cat(str(unique(Results_BP_LM_sel$VAR)))
        cat("\n")
      }
      
      
      ### BOTH INTERVAL & BP exist ---------------------------
      
      LIST_INTERVAL<-list()
      
      ALL_TOGETHER_ENSG_array<-NULL
      
      if(dim(Results_INTERVAL_LM_sel)[1] >0 & dim(Results_BP_LM_sel)[1] >0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("HELLO_WORLD_1\n")
        }
        
        INTERVAL_ENSG_array<-unique(c(Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$CIS_gene_minuslogpvalue)],
                             Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$Block_PCHiC_minuslogpvalue)]))
        
        
        if(Condition_DEBUG == 1)
        {
          cat("INTERVAL_ENSG_array_\n")
          cat(str(INTERVAL_ENSG_array))
          cat("\n")
        }
        
        INTERVAL_DE_gene_df<-data.frame(matrix(ncol=4,nrow=length(INTERVAL_ENSG_array), 
                                               dimnames=list(NULL, c("VAR", "ensembl_gene_id",
                                                                     "Significance",
                                                                     "RNASeq_source"))),
                                        stringsAsFactors = F)
        
        if(length(INTERVAL_ENSG_array) >0)
        {
          ALL_TOGETHER_ENSG_array<-unique(c(ALL_TOGETHER_ENSG_array,INTERVAL_ENSG_array))
          
          
          Results_INTERVAL_LM_sel_ENSG_sel<-Results_INTERVAL_LM_sel[which(Results_INTERVAL_LM_sel$ensembl_gene_id%in%INTERVAL_ENSG_array),]
          
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
          
                   
          
          INTERVAL_DE_gene_df$VAR<-SELECTED_VARS_sel
          INTERVAL_DE_gene_df$ensembl_gene_id<-INTERVAL_ENSG_array
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_DE_gene_df_0\n")
            cat(str(INTERVAL_DE_gene_df))
            cat("\n")
            cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          if(dim(CIS_gene_INTERVAL)[1] >0)
          {
            Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
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
            
            CIS_gene_INTERVAL_subset<-CIS_gene_INTERVAL[,c(which(colnames(CIS_gene_INTERVAL) == "ensembl_gene_id"),
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
            
            INTERVAL_DE_gene_df<-merge(INTERVAL_DE_gene_df,
                              Bind_CIS_Block_INTERVAL,
                              by=c("ensembl_gene_id"),
                              all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("INTERVAL_DE_gene_df_1\n")
              cat(str(INTERVAL_DE_gene_df))
              cat("\n")
              cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
            # quit(status = 1)
            
          }else{
            
            Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
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
            
            INTERVAL_DE_gene_df<-merge(INTERVAL_DE_gene_df,
                              Block_genes_INTERVAL_subset,
                              by=c("ensembl_gene_id"),
                              all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("INTERVAL_DE_gene_df_1\n")
              cat(str(INTERVAL_DE_gene_df))
              cat("\n")
              cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
          }# dim(CIS_gene_INTERVAL)[1] >0
      }# INTERVAL length(INTERVAL_ENSG_array) >0
        
        ############## BP HERE HERE -------------------------
        
        BP_ENSG_array<-unique(c(Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$CIS_gene_minuslogpvalue)],
                             Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$Block_PCHiC_minuslogpvalue)],
                             Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$CIS_gene_minuslogpvalue)],
                             Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$Block_PCHiC_minuslogpvalue)]))
        
       
        if(Condition_DEBUG == 1)
        {
          cat("BP_ENSG_array_\n")
          cat(str(BP_ENSG_array))
          cat("\n")
        }
        
        BP_DE_gene_df<-data.frame(matrix(ncol=4,nrow=length(BP_ENSG_array), 
                                         dimnames=list(NULL, c("VAR", "ensembl_gene_id",
                                                               "Significance",
                                                               "RNASeq_source"))),
                                  stringsAsFactors = F)
        
        if(length(BP_ENSG_array) >0)
        {
          # cat("BP_ENSG_array_\n")
          # cat(str(BP_ENSG_array))
          # cat("\n")
          
         
          
          ALL_TOGETHER_ENSG_array<-unique(c(ALL_TOGETHER_ENSG_array,BP_ENSG_array))
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_BP_LM_sel_\n")
            cat(str(Results_BP_LM_sel))
            cat("\n")
            cat(str(unique(Results_BP_LM_sel$VAR)))
            cat("\n")
          }
          
          Results_BP_LM_sel_ENSG_sel<-Results_BP_LM_sel[which(Results_BP_LM_sel$ensembl_gene_id%in%BP_ENSG_array),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_BP_LM_sel_ENSG_sel_\n")
            cat(str(Results_BP_LM_sel_ENSG_sel))
            cat("\n")
            cat(str(unique(Results_BP_LM_sel_ENSG_sel$ensembl_gene_id)))
            cat("\n")
          }
          
          CIS_gene_BP<-Results_BP_LM_sel_ENSG_sel[!is.na(Results_BP_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                          !is.na(Results_BP_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
          
          if(Condition_DEBUG == 1)
          {
            cat("CIS_gene_\n")
            cat(str(CIS_gene_BP))
            cat("\n")
            cat(str(unique(CIS_gene_BP$ensembl_gene_id)))
            cat("\n")
          }
          
          Block_genes_BP<-Results_BP_LM_sel_ENSG_sel[is.na(Results_BP_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                             !is.na(Results_BP_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Block_genes_\n")
            cat(str(Block_genes_BP))
            cat("\n")
            cat(str(unique(Block_genes_BP$ensembl_gene_id)))
            cat("\n")
          }
          
          #"ajusted.minuslogpvalue_Genotypes","coefficient_Genotypes",
          
          
          
          BP_DE_gene_df$VAR<-SELECTED_VARS_sel
          BP_DE_gene_df$ensembl_gene_id<-BP_ENSG_array
          
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_DE_gene_df_0\n")
            cat(str(BP_DE_gene_df))
            cat("\n")
            cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          if(dim(CIS_gene_BP)[1] >0)
          {
            Block_genes_BP_subset<-Block_genes_BP[,c(which(colnames(Block_genes_BP) == "ensembl_gene_id"),
                                                                 which(colnames(Block_genes_BP) == "HGNC"),
                                                                 which(colnames(Block_genes_BP) == "Cell_Type"),
                                                                 which(colnames(Block_genes_BP) == "Block_PCHiC_minuslogpvalue"),
                                                                 which(colnames(Block_genes_BP) == "Block_PCHiC_Beta"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta_Z_score"))]
            
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
            
            
            if(Condition_DEBUG == 1)
            {
              cat("Block_genes_BP_subset_0\n")
              cat(str(Block_genes_BP_subset))
              cat("\n")
              
            }
            
            CIS_gene_BP_subset<-CIS_gene_BP[,c(which(colnames(CIS_gene_BP) == "ensembl_gene_id"),
                                                           which(colnames(CIS_gene_BP) == "HGNC"),
                                                           which(colnames(CIS_gene_BP) == "Cell_Type"),
                                                           which(colnames(CIS_gene_BP) == "CIS_gene_minuslogpvalue"),
                                                           which(colnames(CIS_gene_BP) == "CIS_Beta"),
                                               which(colnames(CIS_gene_BP) == "CIS_Beta_Z_score"))]
            
            
            colnames(CIS_gene_BP_subset)[which(colnames(CIS_gene_BP_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
            colnames(CIS_gene_BP_subset)[which(colnames(CIS_gene_BP_subset)== "CIS_Beta")]<-"coefficient_Genotypes"
            colnames(CIS_gene_BP_subset)[which(colnames(CIS_gene_BP_subset)== "CIS_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
            
            
            if(Condition_DEBUG == 1)
            {
              cat("CIS_gene_BP_subset_0\n")
              cat(str(CIS_gene_BP_subset))
              cat("\n")
              
            }
            
            
            Bind_CIS_Block_BP<-rbind(CIS_gene_BP_subset,
                                           Block_genes_BP_subset)
            
            if(Condition_DEBUG == 1)
            {
              cat("Bind_CIS_Block_BP_0\n")
              cat(str(Bind_CIS_Block_BP))
              cat("\n")
              
            }
            
            BP_DE_gene_df<-merge(BP_DE_gene_df,
                              Bind_CIS_Block_BP,
                              by=c("ensembl_gene_id"),
                              all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_DE_gene_df_1\n")
              cat(str(BP_DE_gene_df))
              cat("\n")
              cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
            # quit(status = 1)
            
          }else{
            
            Block_genes_BP_subset<-Block_genes_BP[,c(which(colnames(Block_genes_BP) == "ensembl_gene_id"),
                                                     which(colnames(Block_genes_BP) == "HGNC"),
                                                     which(colnames(Block_genes_BP) == "Cell_Type"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_minuslogpvalue"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta_Z_score"))]
            
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
            
            if(Condition_DEBUG == 1)
            {
              cat("Block_genes_BP_subset_0\n")
              cat(str(Block_genes_BP_subset))
              cat("\n")
              
            }
            
            BP_DE_gene_df<-merge(BP_DE_gene_df,
                              Block_genes_BP_subset,
                              by=c("ensembl_gene_id"),
                              all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_DE_gene_df_1\n")
              cat(str(BP_DE_gene_df))
              cat("\n")
              cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
          }# dim(CIS_gene_BP)[1] >0
        }#length(BP_ENSG_array) >0
        
        
        if(length(INTERVAL_ENSG_array) >0 & length(BP_ENSG_array) >0)
        {
          ###################### rbind things together & source ------
          
          
          INTERVAL_DE_gene_df$RNASeq_source<-"Whole blood"
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_DE_gene_df_PRE\n")
            cat(str(INTERVAL_DE_gene_df))
            cat("\n")
            cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          
          
          BP_DE_gene_df$RNASeq_source[which(BP_DE_gene_df$Cell_Type == "Monocyte")]<-"Monocytes"
          BP_DE_gene_df$RNASeq_source[which(BP_DE_gene_df$Cell_Type == "Neutrophil")]<-"Neutrophils"
          BP_DE_gene_df$RNASeq_source[which(BP_DE_gene_df$Cell_Type == "Tcell")]<-"naive T-CD4 Cells"
          
          BP_DE_gene_df<-BP_DE_gene_df[,-which(colnames(BP_DE_gene_df) == "Cell_Type")]
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_DE_gene_df_PRE\n")
            cat(str(BP_DE_gene_df))
            cat("\n")
            cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          
          
          ALL_TOGETHER<-rbind(INTERVAL_DE_gene_df,
                              BP_DE_gene_df)
          
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER_0\n")
            cat(str(ALL_TOGETHER))
            cat("\n")
            cat(str(unique(ALL_TOGETHER$ensembl_gene_id)))
            cat("\n")
          }
          
          ALL_TOGETHER$RNASeq_source<-factor(ALL_TOGETHER$RNASeq_source,
                                             levels=c("Whole blood","Monocytes","Neutrophils","naive T-CD4 Cells"),
                                             ordered=T)
          
          if(Condition_DEBUG == 1)
          {
            cat(sprintf(as.character(names(summary(ALL_TOGETHER$RNASeq_source)))))
            cat("\n")
            cat(sprintf(as.character(summary(ALL_TOGETHER$RNASeq_source))))
            cat("\n")
          }
          
          #### Significance ----
          
          ALL_TOGETHER$Significance[which(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes >= 1.3)]<-"YES"
          ALL_TOGETHER$Significance[which(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes < 1.3)]<-"NO"
          
          ALL_TOGETHER$Significance<-factor(ALL_TOGETHER$Significance,
                                            levels=c("NO","YES"),
                                            ordered=T)
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER$Significance\n")
            cat(sprintf(as.character(names(summary(ALL_TOGETHER$Significance)))))
            cat("\n")
            cat(sprintf(as.character(summary(ALL_TOGETHER$Significance))))
            cat("\n")
          }
          
          
          ### update break.size
          
          local_minuslogpval_max<-max(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes[!is.na(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes)])
          
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
          
          
          
          
          
          ### Use Beta from the full model of beta Z-score ---
          
          ### Labels beta
          
          indx.finite<-is.finite(ALL_TOGETHER$coefficient_Genotypes_Z_score)
          
          check.finite<-sum(indx.finite)
          
          if(check.finite >0)
          {
            A<-summary(ALL_TOGETHER$coefficient_Genotypes_Z_score[indx.finite])
            
            min_value<-A[1]
            
            Q3_value<-A[5]
            max_value<-A[6]
            
            if(Condition_DEBUG == 1)
            {
              cat("Values\n")
              cat(sprintf(as.character((min_value))))
              cat("\n")
              cat(sprintf(as.character((max_value))))
              cat("\n")
              cat(sprintf(as.character((Q3_value))))
              cat("\n")
            }
            
            
            step<-(max_value-min_value)/2
            
            
            if(Condition_DEBUG == 1)
            {
              cat("Step\n")
              cat(sprintf(as.character((step))))
              cat("\n")
            }
            
            candidate_vector<-seq(min_value,max_value+step, by=step)
            
            if(Condition_DEBUG == 1)
            {
              cat("candidate_vector_PRE\n")
              cat(sprintf(as.character((candidate_vector))))
              cat("\n")
            }
            
            if(SELECTED_VARS_sel == "chr17_38764524_T_A")
            {
              # step<-(Q3_value-min_value)/4
              # 
              # candidate_vector<-seq(min_value,Q3_value+step, by=step)
              
              candidate_vector<-seq(-1.75,1.75, by=0.25)
              
              
              if(Condition_DEBUG == 1)
              {
                cat(sprintf(as.character((step))))
                cat("\n")
                cat(sprintf(as.character((candidate_vector))))
                cat("\n")
              }
              
              
            }#SELECTED_VARS_sel == "chr17_38764524_T_A"
            
            if(Condition_DEBUG == 1)
            {
              cat("candidate_vector_DEF\n")
              cat(sprintf(as.character((candidate_vector))))
              cat("\n")
            }
            
            breaks.Beta<-sort(unique(round(c(0,candidate_vector),2)))
            labels.Beta<-as.character(breaks.Beta)
            
          }else{
            
            breaks.Beta<-sort(unique(round(c(0,-2.5,2.5),2)))
            labels.Beta<-as.character(breaks.Beta)
            
          }# length(indx.finite) >0
          
          if(Condition_DEBUG == 1)
          {
            cat("breaks.Beta\n")
            cat(str(breaks.Beta))
            cat("\n")
            
            cat("labels.Beta\n")
            cat(str(labels.Beta))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          ### Graph ----
          
          ALL_TOGETHER<-ALL_TOGETHER[order(ALL_TOGETHER$ensembl_gene_id),]
          
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER_for_dot_plot_2:\t")
            cat(str(ALL_TOGETHER))
            cat("\n")
            
          }
          
          
          
          dotplot<-ggplot(data=ALL_TOGETHER,
                          aes(y=HGNC,
                              x=RNASeq_source)) +
            geom_point(aes(color=Significance,
                           fill=coefficient_Genotypes_Z_score,
                           size=ajusted.minuslogpvalue_Genotypes), stroke=1, shape=21)+
            scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
            scale_size(range = c(0,20), name='-log10pval',
                       breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
            scale_fill_gradient2(
              low = "blue", 
              mid = "white", 
              high = "red", 
              midpoint = 0,
              breaks=breaks.Beta,labels=labels.Beta,
              limits=c(breaks.Beta[1]-0.01,breaks.Beta[length(breaks.Beta)]+0.01),name=paste('Effect size','Z score', sep="\n"),na.value = "gray")+
            scale_y_discrete(name=NULL, drop=F)+
            theme_classic()+
            ggeasy::easy_center_title()
          
          dotplot<-dotplot+
            facet_grid(cols = vars(ALL_TOGETHER$RNASeq_source), scales='free_x', space='free_x',
                       drop=F) +
            theme_cowplot(font_size = 14)+
            theme( strip.background = element_blank(),
                   strip.placement = "inside",
                   strip.text = element_text(size=14, angle=90),
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
          
          
          setwd(path8)
          # setwd(out)
          
          svgname<-paste("DE_ANALYSIS_",SELECTED_VARS_sel,".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= dotplot,
                   device="svg",
                   height=10, width=12)
          }
          
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("PRINTING#1\n")
            
          }
          # ############################################### HERE HERE 
          # quit(status = 1)
          
          
        }#length(INTERVAL_ENSG_array) >0 & length(BP_ENSG_array) >0
       
      }#dim(Results_INTERVAL_LM_sel)[1] >0 & dim(Results_BP_LM_sel)[1] >0
      
      if(dim(Results_INTERVAL_LM_sel)[1] >0 & dim(Results_BP_LM_sel)[1] == 0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("HELLO_WORLD_2\n")
        }
        
        INTERVAL_ENSG_array<-unique(c(Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$CIS_gene_minuslogpvalue)],
                                      Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$Block_PCHiC_minuslogpvalue)]))
        
        if(Condition_DEBUG == 1)
        {
          cat("INTERVAL_ENSG_array_\n")
          cat(str(INTERVAL_ENSG_array))
          cat("\n")
        }
        
        INTERVAL_DE_gene_df<-data.frame(matrix(ncol=4,nrow=length(INTERVAL_ENSG_array), 
                                               dimnames=list(NULL, c("VAR", "ensembl_gene_id",
                                                                     "Significance",
                                                                     "RNASeq_source"))),
                                        stringsAsFactors = F) 
       
       
        if(length(INTERVAL_ENSG_array) >0)
        {
          
          
          
          ALL_TOGETHER_ENSG_array<-unique(c(ALL_TOGETHER_ENSG_array,INTERVAL_ENSG_array))
          
          
          Results_INTERVAL_LM_sel_ENSG_sel<-Results_INTERVAL_LM_sel[which(Results_INTERVAL_LM_sel$ensembl_gene_id%in%INTERVAL_ENSG_array),]
          
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
          
          
          
          INTERVAL_DE_gene_df$VAR<-SELECTED_VARS_sel
          INTERVAL_DE_gene_df$ensembl_gene_id<-INTERVAL_ENSG_array
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_DE_gene_df_0\n")
            cat(str(INTERVAL_DE_gene_df))
            cat("\n")
            cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          if(dim(CIS_gene_INTERVAL)[1] >0)
          {
            Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
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
            
            CIS_gene_INTERVAL_subset<-CIS_gene_INTERVAL[,c(which(colnames(CIS_gene_INTERVAL) == "ensembl_gene_id"),
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
            
            INTERVAL_DE_gene_df<-merge(INTERVAL_DE_gene_df,
                                       Bind_CIS_Block_INTERVAL,
                                       by=c("ensembl_gene_id"),
                                       all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("INTERVAL_DE_gene_df_1\n")
              cat(str(INTERVAL_DE_gene_df))
              cat("\n")
              cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
            # quit(status = 1)
            
          }else{
            
            Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
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
            
            INTERVAL_DE_gene_df<-merge(INTERVAL_DE_gene_df,
                                       Block_genes_INTERVAL_subset,
                                       by=c("ensembl_gene_id"),
                                       all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("INTERVAL_DE_gene_df_1\n")
              cat(str(INTERVAL_DE_gene_df))
              cat("\n")
              cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
          }# dim(CIS_gene_INTERVAL)[1] >0
        }# INTERVAL length(INTERVAL_ENSG_array) >0
        
    
        
        if(length(INTERVAL_ENSG_array) >0)
        {
          ###################### rbind things together & source ------
          
          
          INTERVAL_DE_gene_df$RNASeq_source<-"Whole blood"
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_DE_gene_df_PRE\n")
            cat(str(INTERVAL_DE_gene_df))
            cat("\n")
            cat(str(unique(INTERVAL_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          
          
          
          
          ALL_TOGETHER<-rbind(INTERVAL_DE_gene_df)
          
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER_0\n")
            cat(str(ALL_TOGETHER))
            cat("\n")
            cat(str(unique(ALL_TOGETHER$ensembl_gene_id)))
            cat("\n")
          }
          
          ALL_TOGETHER$RNASeq_source<-factor(ALL_TOGETHER$RNASeq_source,
                                             levels=c("Whole blood","Monocytes","Neutrophils","naive T-CD4 Cells"),
                                             ordered=T)
          
          if(Condition_DEBUG == 1)
          {
            cat(sprintf(as.character(names(summary(ALL_TOGETHER$RNASeq_source)))))
            cat("\n")
            cat(sprintf(as.character(summary(ALL_TOGETHER$RNASeq_source))))
            cat("\n")
          }
          
          #### Significance ----
          
          ALL_TOGETHER$Significance[which(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes >= 1.3)]<-"YES"
          ALL_TOGETHER$Significance[which(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes < 1.3)]<-"NO"
          
          ALL_TOGETHER$Significance<-factor(ALL_TOGETHER$Significance,
                                            levels=c("NO","YES"),
                                            ordered=T)
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER$Significance\n")
            cat(sprintf(as.character(names(summary(ALL_TOGETHER$Significance)))))
            cat("\n")
            cat(sprintf(as.character(summary(ALL_TOGETHER$Significance))))
            cat("\n")
          }
          
          
          ### update break.size
          
          local_minuslogpval_max<-max(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes[!is.na(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes)])
          
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
          
          
          
          
          
          ### Use Beta from the full model of beta Z-score ---
          
          ### Labels beta
          
          indx.finite<-is.finite(ALL_TOGETHER$coefficient_Genotypes_Z_score)
          
          check.finite<-sum(indx.finite)
          
          if(check.finite >0)
          {
            A<-summary(ALL_TOGETHER$coefficient_Genotypes_Z_score[indx.finite])
            
            min_value<-A[1]
            
            Q3_value<-A[5]
            max_value<-A[6]
            
            if(Condition_DEBUG == 1)
            {
              cat("Values\n")
              cat(sprintf(as.character((min_value))))
              cat("\n")
              cat(sprintf(as.character((max_value))))
              cat("\n")
              cat(sprintf(as.character((Q3_value))))
              cat("\n")
            }
            
            
            step<-(max_value-min_value)/2
            
            
            if(Condition_DEBUG == 1)
            {
              cat("Step\n")
              cat(sprintf(as.character((step))))
              cat("\n")
            }
            
            candidate_vector<-seq(min_value,max_value+step, by=step)
            
            if(Condition_DEBUG == 1)
            {
              cat("candidate_vector_PRE\n")
              cat(sprintf(as.character((candidate_vector))))
              cat("\n")
            }
            
            
            
            
            breaks.Beta<-sort(unique(round(c(0,candidate_vector),2)))
            labels.Beta<-as.character(breaks.Beta)
            
          }else{
            
            breaks.Beta<-sort(unique(round(c(0,-2.5,2.5),2)))
            labels.Beta<-as.character(breaks.Beta)
            
          }# length(indx.finite) >0
          
          if(Condition_DEBUG == 1)
          {
            cat("breaks.Beta\n")
            cat(str(breaks.Beta))
            cat("\n")
            
            cat("labels.Beta\n")
            cat(str(labels.Beta))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          ### Graph ----
          
          ALL_TOGETHER<-ALL_TOGETHER[order(ALL_TOGETHER$ensembl_gene_id),]
          
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER_for_dot_plot_2:\t")
            cat(str(ALL_TOGETHER))
            cat("\n")
            
          }
          
          
          
          dotplot<-ggplot(data=ALL_TOGETHER,
                          aes(y=HGNC,
                              x=RNASeq_source)) +
            geom_point(aes(color=Significance,
                           fill=coefficient_Genotypes_Z_score,
                           size=ajusted.minuslogpvalue_Genotypes), stroke=1, shape=21)+
            scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
            scale_size(range = c(0,20), name='-log10pval',
                       breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
            scale_fill_gradient2(
              low = "blue", 
              mid = "white", 
              high = "red", 
              midpoint = 0,
              breaks=breaks.Beta,labels=labels.Beta,
              limits=c(breaks.Beta[1]-0.01,breaks.Beta[length(breaks.Beta)]+0.01),name=paste('Effect size','Z score', sep="\n"),na.value = "gray")+
            scale_y_discrete(name=NULL, drop=F)+
            theme_classic()+
            scale_x_discrete(name=NULL, drop=T)+
            ggeasy::easy_center_title()
          
          dotplot<-dotplot+
            facet_grid(cols = vars(ALL_TOGETHER$RNASeq_source), scales='free_x', space='free_x',
                       drop=T) +
            theme_cowplot(font_size = 14)+
            theme( strip.background = element_blank(),
                   strip.placement = "inside",
                   strip.text = element_text(size=14, angle=90),
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
                  legend.text = element_text(size=14)) #change legend text font size
           
          
          setwd(path8)
          # setwd(out)
          
          svgname<-paste("DE_ANALYSIS_",SELECTED_VARS_sel,".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= dotplot,
                   device="svg",
                   height=10, width=12)
          }
          
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("PRINTING#CASE2#\n")
            
          }
          # ############################################### HERE HERE 
          # quit(status = 1)
          
        }#length(INTERVAL_ENSG_array) >0
      }#dim(Results_INTERVAL_LM_sel)[1] >0 & dim(Results_BP_LM_sel)[1] == 0
      
      if(dim(Results_INTERVAL_LM_sel)[1] ==0 & dim(Results_BP_LM_sel)[1] >0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("HELLO_WORLD_3\n")
        }
               
        BP_ENSG_array<-unique(c(Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$CIS_gene_minuslogpvalue)],
                                Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$Block_PCHiC_minuslogpvalue)],
                                Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$CIS_gene_minuslogpvalue)],
                                Results_BP_LM_sel$ensembl_gene_id[!is.na(Results_BP_LM_sel$Block_PCHiC_minuslogpvalue)]))
        
        if(Condition_DEBUG == 1)
        {
          cat("BP_ENSG_array_\n")
          cat(str(BP_ENSG_array))
          cat("\n")
        }
        
        BP_DE_gene_df<-data.frame(matrix(ncol=4,nrow=length(BP_ENSG_array), 
                                         dimnames=list(NULL, c("VAR", "ensembl_gene_id",
                                                               "Significance",
                                                               "RNASeq_source"))),
                                  stringsAsFactors = F)
        
        if(length(BP_ENSG_array) >0)
        {
         
          
        
          ALL_TOGETHER_ENSG_array<-unique(c(ALL_TOGETHER_ENSG_array,BP_ENSG_array))
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_BP_LM_sel_\n")
            cat(str(Results_BP_LM_sel))
            cat("\n")
            cat(str(unique(Results_BP_LM_sel$VAR)))
            cat("\n")
          }
          
          Results_BP_LM_sel_ENSG_sel<-Results_BP_LM_sel[which(Results_BP_LM_sel$ensembl_gene_id%in%BP_ENSG_array),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_BP_LM_sel_ENSG_sel_\n")
            cat(str(Results_BP_LM_sel_ENSG_sel))
            cat("\n")
            cat(str(unique(Results_BP_LM_sel_ENSG_sel$ensembl_gene_id)))
            cat("\n")
          }
          
          CIS_gene_BP<-Results_BP_LM_sel_ENSG_sel[!is.na(Results_BP_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                    !is.na(Results_BP_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
          
          if(Condition_DEBUG == 1)
          {
            cat("CIS_gene_\n")
            cat(str(CIS_gene_BP))
            cat("\n")
            cat(str(unique(CIS_gene_BP$ensembl_gene_id)))
            cat("\n")
          }
          
          Block_genes_BP<-Results_BP_LM_sel_ENSG_sel[is.na(Results_BP_LM_sel_ENSG_sel$CIS_gene_minuslogpvalue) &
                                                       !is.na(Results_BP_LM_sel_ENSG_sel$Block_PCHiC_minuslogpvalue),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Block_genes_\n")
            cat(str(Block_genes_BP))
            cat("\n")
            cat(str(unique(Block_genes_BP$ensembl_gene_id)))
            cat("\n")
          }
          
          #"ajusted.minuslogpvalue_Genotypes","coefficient_Genotypes",
          
          
          
          BP_DE_gene_df$VAR<-SELECTED_VARS_sel
          BP_DE_gene_df$ensembl_gene_id<-BP_ENSG_array
          
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_DE_gene_df_0\n")
            cat(str(BP_DE_gene_df))
            cat("\n")
            cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          if(dim(CIS_gene_BP)[1] >0)
          {
            Block_genes_BP_subset<-Block_genes_BP[,c(which(colnames(Block_genes_BP) == "ensembl_gene_id"),
                                                     which(colnames(Block_genes_BP) == "HGNC"),
                                                     which(colnames(Block_genes_BP) == "Cell_Type"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_minuslogpvalue"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta_Z_score"))]
            
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
            
            
            if(Condition_DEBUG == 1)
            {
              cat("Block_genes_BP_subset_0\n")
              cat(str(Block_genes_BP_subset))
              cat("\n")
              
            }
            
            CIS_gene_BP_subset<-CIS_gene_BP[,c(which(colnames(CIS_gene_BP) == "ensembl_gene_id"),
                                               which(colnames(CIS_gene_BP) == "HGNC"),
                                               which(colnames(CIS_gene_BP) == "Cell_Type"),
                                               which(colnames(CIS_gene_BP) == "CIS_gene_minuslogpvalue"),
                                               which(colnames(CIS_gene_BP) == "CIS_Beta"),
                                               which(colnames(CIS_gene_BP) == "CIS_Beta_Z_score"))]
            
            
            colnames(CIS_gene_BP_subset)[which(colnames(CIS_gene_BP_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
            colnames(CIS_gene_BP_subset)[which(colnames(CIS_gene_BP_subset)== "CIS_Beta")]<-"coefficient_Genotypes"
            colnames(CIS_gene_BP_subset)[which(colnames(CIS_gene_BP_subset)== "CIS_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
            
            
            if(Condition_DEBUG == 1)
            {
              cat("CIS_gene_BP_subset_0\n")
              cat(str(CIS_gene_BP_subset))
              cat("\n")
              
            }
            
            
            Bind_CIS_Block_BP<-rbind(CIS_gene_BP_subset,
                                     Block_genes_BP_subset)
            
            if(Condition_DEBUG == 1)
            {
              cat("Bind_CIS_Block_BP_0\n")
              cat(str(Bind_CIS_Block_BP))
              cat("\n")
              
            }
            
            BP_DE_gene_df<-merge(BP_DE_gene_df,
                                 Bind_CIS_Block_BP,
                                 by=c("ensembl_gene_id"),
                                 all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_DE_gene_df_1\n")
              cat(str(BP_DE_gene_df))
              cat("\n")
              cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
            # quit(status = 1)
            
          }else{
            
            Block_genes_BP_subset<-Block_genes_BP[,c(which(colnames(Block_genes_BP) == "ensembl_gene_id"),
                                                     which(colnames(Block_genes_BP) == "HGNC"),
                                                     which(colnames(Block_genes_BP) == "Cell_Type"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_minuslogpvalue"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta"),
                                                     which(colnames(Block_genes_BP) == "Block_PCHiC_Beta_Z_score"))]
            
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta")]<-"coefficient_Genotypes"
            colnames(Block_genes_BP_subset)[which(colnames(Block_genes_BP_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Genotypes_Z_score"
            
            if(Condition_DEBUG == 1)
            {
              cat("Block_genes_BP_subset_0\n")
              cat(str(Block_genes_BP_subset))
              cat("\n")
              
            }
            
            BP_DE_gene_df<-merge(BP_DE_gene_df,
                                 Block_genes_BP_subset,
                                 by=c("ensembl_gene_id"),
                                 all=T)
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_DE_gene_df_1\n")
              cat(str(BP_DE_gene_df))
              cat("\n")
              cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
              cat("\n")
            }
            
          }# dim(CIS_gene_BP)[1] >0
        }#length(BP_ENSG_array) >0
        
        if(length(BP_ENSG_array) >0)
        {
          ###################### rbind things together & source ------
          
          
          
          
          
          BP_DE_gene_df$RNASeq_source[which(BP_DE_gene_df$Cell_Type == "Monocyte")]<-"Monocytes"
          BP_DE_gene_df$RNASeq_source[which(BP_DE_gene_df$Cell_Type == "Neutrophil")]<-"Neutrophils"
          BP_DE_gene_df$RNASeq_source[which(BP_DE_gene_df$Cell_Type == "Tcell")]<-"naive T-CD4 Cells"
          
          BP_DE_gene_df<-BP_DE_gene_df[,-which(colnames(BP_DE_gene_df) == "Cell_Type")]
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_DE_gene_df_PRE\n")
            cat(str(BP_DE_gene_df))
            cat("\n")
            cat(str(unique(BP_DE_gene_df$ensembl_gene_id)))
            cat("\n")
          }
          
          
          
          ALL_TOGETHER<-rbind(BP_DE_gene_df)
          
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER_0\n")
            cat(str(ALL_TOGETHER))
            cat("\n")
            cat(str(unique(ALL_TOGETHER$ensembl_gene_id)))
            cat("\n")
          }
          
          ALL_TOGETHER$RNASeq_source<-factor(ALL_TOGETHER$RNASeq_source,
                                             levels=c("Whole blood","Monocytes","Neutrophils","naive T-CD4 Cells"),
                                             ordered=T)
          
          if(Condition_DEBUG == 1)
          {
            cat(sprintf(as.character(names(summary(ALL_TOGETHER$RNASeq_source)))))
            cat("\n")
            cat(sprintf(as.character(summary(ALL_TOGETHER$RNASeq_source))))
            cat("\n")
          }
          
          #### Significance ----
          
          ALL_TOGETHER$Significance[which(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes >= 1.3)]<-"YES"
          ALL_TOGETHER$Significance[which(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes < 1.3)]<-"NO"
          
          ALL_TOGETHER$Significance<-factor(ALL_TOGETHER$Significance,
                                            levels=c("NO","YES"),
                                            ordered=T)
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER$Significance\n")
            cat(sprintf(as.character(names(summary(ALL_TOGETHER$Significance)))))
            cat("\n")
            cat(sprintf(as.character(summary(ALL_TOGETHER$Significance))))
            cat("\n")
          }
          
          
          ### update break.size
          
          local_minuslogpval_max<-max(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes[!is.na(ALL_TOGETHER$ajusted.minuslogpvalue_Genotypes)])
          
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
          
          
          
          
          
          ### Use Beta from the full model of beta Z-score ---
          
          ### Labels beta
          
          indx.finite<-is.finite(ALL_TOGETHER$coefficient_Genotypes_Z_score)
          
          check.finite<-sum(indx.finite)
          
          if(check.finite >0)
          {
            A<-summary(ALL_TOGETHER$coefficient_Genotypes_Z_score[indx.finite])
            
            min_value<-A[1]
            
            Q3_value<-A[5]
            max_value<-A[6]
            
            if(Condition_DEBUG == 1)
            {
              cat("Values\n")
              cat(sprintf(as.character((min_value))))
              cat("\n")
              cat(sprintf(as.character((max_value))))
              cat("\n")
              cat(sprintf(as.character((Q3_value))))
              cat("\n")
            }
            
            
            step<-(max_value-min_value)/2
            
            
            if(Condition_DEBUG == 1)
            {
              cat("Step\n")
              cat(sprintf(as.character((step))))
              cat("\n")
            }
            
            candidate_vector<-seq(min_value,max_value+step, by=step)
            
            if(Condition_DEBUG == 1)
            {
              cat("candidate_vector_PRE\n")
              cat(sprintf(as.character((candidate_vector))))
              cat("\n")
            }
            
            if(SELECTED_VARS_sel == "chr17_38764524_T_A")
            {
              # step<-(Q3_value-min_value)/4
              # 
              # candidate_vector<-seq(min_value,Q3_value+step, by=step)
              
              candidate_vector<-seq(-1.75,1.75, by=0.25)
              
              
              if(Condition_DEBUG == 1)
              {
                cat(sprintf(as.character((step))))
                cat("\n")
                cat(sprintf(as.character((candidate_vector))))
                cat("\n")
              }
              
              
            }#SELECTED_VARS_sel == "chr17_38764524_T_A"
            
            if(Condition_DEBUG == 1)
            {
              cat("candidate_vector_DEF\n")
              cat(sprintf(as.character((candidate_vector))))
              cat("\n")
            }
            
            breaks.Beta<-sort(unique(round(c(0,candidate_vector),2)))
            labels.Beta<-as.character(breaks.Beta)
            
          }else{
            
            breaks.Beta<-sort(unique(round(c(0,-2.5,2.5),2)))
            labels.Beta<-as.character(breaks.Beta)
            
          }# length(indx.finite) >0
          
          if(Condition_DEBUG == 1)
          {
            cat("breaks.Beta\n")
            cat(str(breaks.Beta))
            cat("\n")
            
            cat("labels.Beta\n")
            cat(str(labels.Beta))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          ### Graph ----
          
          ALL_TOGETHER<-ALL_TOGETHER[order(ALL_TOGETHER$ensembl_gene_id),]
          
          if(Condition_DEBUG == 1)
          {
            cat("ALL_TOGETHER_for_dot_plot_2:\t")
            cat(str(ALL_TOGETHER))
            cat("\n")
            
          }
          
          
          
          dotplot<-ggplot(data=ALL_TOGETHER,
                          aes(y=HGNC,
                              x=RNASeq_source)) +
            geom_point(aes(color=Significance,
                           fill=coefficient_Genotypes_Z_score,
                           size=ajusted.minuslogpvalue_Genotypes), stroke=1, shape=21)+
            scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
            scale_size(range = c(0,20), name='-log10pval',
                       breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
            scale_fill_gradient2(
              low = "blue", 
              mid = "white", 
              high = "red", 
              midpoint = 0,
              breaks=breaks.Beta,labels=labels.Beta,
              limits=c(breaks.Beta[1]-0.01,breaks.Beta[length(breaks.Beta)]+0.01),name=paste('Effect size','Z score', sep="\n"),na.value = "gray")+
            scale_y_discrete(name=NULL, drop=F)+
            theme_classic()+
            scale_x_discrete(name=NULL, drop=T)+
            ggeasy::easy_center_title()
          
          dotplot<-dotplot+
            facet_grid(cols = vars(ALL_TOGETHER$RNASeq_source), scales='free_x', space='free_x',
                       drop=T) +
            theme_cowplot(font_size = 14)+
            theme( strip.background = element_blank(),
                   strip.placement = "inside",
                   strip.text = element_text(size=14, angle=90),
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
            
          
          
          setwd(path8)
          # setwd(out)
          
          svgname<-paste("DE_ANALYSIS_",SELECTED_VARS_sel,".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= dotplot,
                   device="svg",
                   height=10, width=12)
          }
          
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("PRINTING#1\n")
            
          }
          # ############################################### HERE HERE 
          # quit(status = 1)
          
        }#length(BP_ENSG_array) >0
        
       
      }#dim(Results_INTERVAL_LM_sel)[1] ==0 & dim(Results_BP_LM_sel)[1] >0
      ### Neither INTERVAL nor BP exist ----
      
      if(dim(Results_INTERVAL_LM_sel)[1] == 0 & dim(Results_BP_LM_sel)[1] == 0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("HELLO_WORLD_4\n")
        }
        
        # quit(status = 1)
      }
      
      
      # ################################## HERE HERE HERE  GET THE drop=T back to drop=F and comment setwd(out)
      # quit(status = 1)
      
      
      ############################# VIOLIN PLOTS ---------------------------------------------------------------------------------------------------------------------------------------
      if(Condition_DEBUG == 1)
      {
        cat("ALL_TOGETHER_ENSG_array_FOR_VIOLIN_PLOTS\n")
        cat(str(ALL_TOGETHER_ENSG_array))
        cat("\n")
      }
      
      if(length(ALL_TOGETHER_ENSG_array) >0)
      {
        Condition_DEBUG <- 0
        
        for(z in 1:length(ALL_TOGETHER_ENSG_array))
        {
          
          ALL_TOGETHER_ENSG_array_sel<-ALL_TOGETHER_ENSG_array[z]
          
          cat("---------------->\t")
          cat(sprintf(as.character(ALL_TOGETHER_ENSG_array_sel)))
          cat("\t")
          
          ALL_TOGETHER_ENSG_sel<-ALL_TOGETHER[which(ALL_TOGETHER$ensembl_gene_id%in%ALL_TOGETHER_ENSG_array_sel),]
          
          HGNC_sel<-unique(ALL_TOGETHER_ENSG_sel$HGNC)
          
          cat(sprintf(as.character(HGNC_sel)))
          cat("\n")
          
          if(Condition_DEBUG == 1)
          {
            cat("residuals_INTERVAL_df\n")
            cat(str(residuals_INTERVAL_df))
            cat("\n")
            
          }
          
          residuals_INTERVAL_df_sel_ENSG_sel<-residuals_INTERVAL_df[which(residuals_INTERVAL_df$ensembl_gene_id == ALL_TOGETHER_ENSG_array_sel),]
          
          if(Condition_DEBUG == 1)
          {
            cat("residuals_INTERVAL_df_sel_ENSG_sel_\n")
            cat(str(residuals_INTERVAL_df_sel_ENSG_sel))
            cat("\n")
            
          }
          
          INTERVAL_covariates_and_PEER_factors_sel_subset<-INTERVAL_covariates_and_PEER_factors_sel[,c(which(colnames(INTERVAL_covariates_and_PEER_factors_sel) == "Genotype"),
                                                                                                       which(colnames(INTERVAL_covariates_and_PEER_factors_sel) == "sample_id"))]
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_covariates_and_PEER_factors_sel_subset_\n")
            cat(str(INTERVAL_covariates_and_PEER_factors_sel_subset))
            cat("\n")
            
          }
          
          INTERVAL_GENE_EXP_sel<-INTERVAL_GENE_EXP[,c(which(colnames(INTERVAL_GENE_EXP) == "sample_id"),
                                                      which(colnames(INTERVAL_GENE_EXP) %in% ALL_TOGETHER_ENSG_array_sel))]
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_GENE_EXP_sel_1\n")
            cat(str(INTERVAL_GENE_EXP_sel))
            cat("\n")
          }
          
          residuals_INTERVAL_df_sel_ENSG_sel.m<-melt(residuals_INTERVAL_df_sel_ENSG_sel, id.vars=c("ensembl_gene_id"),value.name = "residuals_FPKM",variable.name = "sample_id")
          
          
          if(Condition_DEBUG == 1)
          {
            # cat("residuals_INTERVAL_df_sel_ENSG_sel\n")
            # cat(str(residuals_INTERVAL_df_sel_ENSG_sel))
            # cat("\n")
            
            cat("residuals_INTERVAL_df_sel_ENSG_sel.m\n")
            cat(str(residuals_INTERVAL_df_sel_ENSG_sel.m))
            cat("\n")
          }
          
          ################# melt the INTERVAL_GENE_EXP_sel
          
          INTERVAL_GENE_EXP_sel.m<-melt(INTERVAL_GENE_EXP_sel, id.vars=c("sample_id"),value.name = "log2FPKM",variable.name = "ensembl_gene_id")
          INTERVAL_GENE_EXP_sel.m$FPKM<-2^(INTERVAL_GENE_EXP_sel.m$log2FPKM)
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_GENE_EXP_sel.m_0\n")
            cat(str(INTERVAL_GENE_EXP_sel.m))
            cat("\n")
            
            cat("sample_id\n")
            cat(str(INTERVAL_GENE_EXP_sel.m$sample_id))
            cat("\n")
            
            cat("ensembl_gene_id\n")
            cat(str(INTERVAL_GENE_EXP_sel.m$ensembl_gene_id))
            cat("\n")
            
            # ########################################
            # quit(status = 1)
          }
          
          INTERVAL_GENE_EXP_sel.m<-merge(INTERVAL_GENE_EXP_sel.m,
                                         residuals_INTERVAL_df_sel_ENSG_sel.m,
                                         by=c("sample_id","ensembl_gene_id"))
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_GENE_EXP_sel.m_1\n")
            cat(str(INTERVAL_GENE_EXP_sel.m))
            cat("\n")
            
            # #####################################################
            # quit(status = 1)
          }
          
          A<-round(summary(INTERVAL_GENE_EXP_sel.m$residuals_FPKM[!is.na(INTERVAL_GENE_EXP_sel.m$residuals_FPKM)]),2)
          # A2<-round(summary(INTERVAL_GENE_EXP_sel.m$FPKM[!is.na(INTERVAL_GENE_EXP_sel.m$FPKM)]),2)
          
          INTERVAL_GENE_EXP_sel.m$residuals_FPKM_no_negative<-INTERVAL_GENE_EXP_sel.m$residuals_FPKM+abs(A[1])
          
          A3<-round(summary(INTERVAL_GENE_EXP_sel.m$residuals_FPKM_no_negative[!is.na(INTERVAL_GENE_EXP_sel.m$residuals_FPKM_no_negative)]),0)
          
          
          
          if(Condition_DEBUG == 1)
          {
            cat("INTERVAL_GENE_EXP_sel.m_2\n")
            cat(str(INTERVAL_GENE_EXP_sel.m))
            cat("\n")
            
            cat("Summary_residuals_FPKM:\t")
            cat(sprintf(as.character(names(A))))
            cat("\n")
            cat(sprintf(as.character(A)))
            cat("\n")
            
            # cat("Summary_FPKM:\t")
            # cat(sprintf(as.character(names(A2))))
            # cat("\n")
            # cat(sprintf(as.character(A2)))
            # cat("\n")
            
            
            cat("Summary_residuals_FPKM_no_negative:\t")
            cat(sprintf(as.character(names(A3))))
            cat("\n")
            cat(sprintf(as.character(A3)))
            cat("\n")
            
            
          }
          
          
          
          #### Merge with covariates and Keep to HET ----
          
          
          
          
          INTERVAL_GENE_EXP_sel.m<-merge(INTERVAL_GENE_EXP_sel.m,
                                         INTERVAL_covariates_and_PEER_factors_sel_subset,
                                         by=c("sample_id"))
          
          
          if(Condition_DEBUG == 1)
          {
            
            cat("INTERVAL_GENE_EXP_sel.m_3\n")
            cat(str(INTERVAL_GENE_EXP_sel.m))
            cat("\n")
            
            # quit(status = 1)
            
          }
          
          INTERVAL_GENE_EXP_sel.m_HET<-droplevels(INTERVAL_GENE_EXP_sel.m[which(INTERVAL_GENE_EXP_sel.m$Genotype != "HOM"),])
          
          if(Condition_DEBUG == 1)
          {
            
            cat("INTERVAL_GENE_EXP_sel.m_HET_0\n")
            cat(str(INTERVAL_GENE_EXP_sel.m_HET))
            cat("\n")
            
            
            # quit(status = 1)
            
          }
          
         
          
          ##### Violin visual log2FPKM unadjusted----
          
          A<-round(summary(INTERVAL_GENE_EXP_sel.m_HET$log2FPKM[!is.na(INTERVAL_GENE_EXP_sel.m_HET$log2FPKM)]),2)
          step<-abs(A[6]-A[1])/10
          
          if(Condition_DEBUG == 1)
          {
            cat("summary_INTERVAL_GENE_EXP_sel.m_HET$log2FPK\n")
            cat(sprintf(as.character(names(A))))
            cat("\n")
            cat(sprintf(as.character(A)))
            cat("\n")
          }
          
          if(step == 0)
          {
            
            step<-1
          }
          
          breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
          labels.Rank<-as.character(round(breaks.Rank,3))
          
          levels_ensembl_gene_id<-unique(INTERVAL_GENE_EXP_sel.m_HET$ensembl_gene_id)
          
          if(Condition_DEBUG == 1)
          {
            cat("labels.Rank:\t")
            cat(sprintf(as.character(labels.Rank)))
            cat("\n")
            
            cat("levels_ensembl_gene_id:\t")
            cat(sprintf(as.character(levels_ensembl_gene_id)))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          
          
          graph_log2FPKM<-ggplot(data=INTERVAL_GENE_EXP_sel.m_HET,aes(x=Genotype, y=log2FPKM, fill=Genotype)) +
            geom_violin()+
            stat_summary(fun = median, fun.min = median, fun.max = median,
                         geom = "crossbar", width = 0.5)+
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_blank(),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_discrete(name=NULL, drop=F)+
            scale_y_continuous(name="Gene EXP normalised (log2FPKM)",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
            scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=12))+
            ggeasy::easy_center_title()
          
          
          ##### Violin visual FPKM_adjusted reduced model ----
          
          INTERVAL_GENE_EXP_sel.m_HET$Log2_residuals_FPKM_no_negative<-log2(INTERVAL_GENE_EXP_sel.m_HET$residuals_FPKM_no_negative)
          
          if(Condition_DEBUG == 1)
          {
            
            cat("INTERVAL_GENE_EXP_sel.m_HET_FPKM_reduced_model\n")
            cat(str(INTERVAL_GENE_EXP_sel.m_HET))
            cat("\n")
            
            
            # quit(status = 1)
            
          }
          
          
          A<-round(summary(INTERVAL_GENE_EXP_sel.m_HET$Log2_residuals_FPKM_no_negative[!is.na(INTERVAL_GENE_EXP_sel.m_HET$Log2_residuals_FPKM_no_negative)]),2)
         
          step<-abs(A[6]-A[1])/10
          
          if(Condition_DEBUG == 1)
          {
            cat("summary_INTERVAL_GENE_EXP_sel.m_HET$Log2_residuals_FPKM_no_negative\n")
            cat(sprintf(as.character(names(A))))
            cat("\n")
            cat(sprintf(as.character(A)))
            cat("\n")
          }
          
          if(step == 0)
          {
            
            step<-1
          }
          
          breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
          labels.Rank<-as.character(round(breaks.Rank,3))
          
          levels_ensembl_gene_id<-unique(INTERVAL_GENE_EXP_sel.m_HET$ensembl_gene_id)
          
          if(Condition_DEBUG == 1)
          {
            cat("labels.Rank:\t")
            cat(sprintf(as.character(labels.Rank)))
            cat("\n")
            
            cat("levels_ensembl_gene_id:\t")
            cat(sprintf(as.character(levels_ensembl_gene_id)))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          
          
          graph_adjusted_FPKM<-ggplot(data=INTERVAL_GENE_EXP_sel.m_HET,aes(x=Genotype, y=Log2_residuals_FPKM_no_negative, fill=Genotype)) +
            geom_violin()+
            stat_summary(fun = median, fun.min = median, fun.max = median,
                         geom = "crossbar", width = 0.5)+
            theme_bw()+
            theme(axis.title.y=element_text(size=24, family="sans"),
                  axis.title.x=element_text(size=24, family="sans"),
                  axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                  axis.text.x=element_blank(),
                  legend.title=element_text(size=16,color="black", family="sans"),
                  legend.text=element_text(size=12,color="black", family="sans"))+
            scale_x_discrete(name=NULL, drop=F)+
            scale_y_continuous(name="Residuals model Gene EXP log2FPKM adjusted for covariates",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
            scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
            theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=12))+
            ggeasy::easy_center_title()
          
          ################################# BP EXP if any --------------------------
          
          path_BP_VAR<-paste(path_BP,SELECTED_VARS_sel,'/',sep='')
          
          if(Condition_DEBUG == 1)
          {
            cat("path_BP_VAR:\t")
            cat(sprintf(as.character(path_BP_VAR)))
            cat("\n")
            
            # quit(status = 1)
          }
          
          setwd(path_BP_VAR)
          
          FLAG_BP<-"NA"
          
          graph_log2FPKM_BP<-NULL
          
          if (file.exists("BP_for_LM.rds")){
            
            FLAG_BP<-1
            BP_genes_DEF<-readRDS(file="BP_for_LM.rds")
            
            if(Condition_DEBUG == 1)
            {
              
              cat("BP_genes_DEF_0\n")
              cat(str(BP_genes_DEF))
              cat("\n")
            }
            
            BP_genes_DEF_sel.m<-unique(BP_genes_DEF[which(BP_genes_DEF$ensembl_gene_id == ALL_TOGETHER_ENSG_array_sel),c(which(colnames(BP_genes_DEF) == "ensembl_gene_id"),
                                                                                                                         which(colnames(BP_genes_DEF) == "HGNC"),
                                                                                                                         which(colnames(BP_genes_DEF) == "value"),
                                                                                                                         which(colnames(BP_genes_DEF) == "Cell_Type"),
                                                                                                                         which(colnames(BP_genes_DEF) == "Genotype"))])
            if(dim(BP_genes_DEF_sel.m)[1] >0)
            {
              
              HGNC_sel<-unique(BP_genes_DEF_sel.m$HGNC)
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_genes_DEF_sel.m_0\n")
                cat(str(BP_genes_DEF_sel.m))
                cat("\n")
              }
              
              
              
              BP_genes_DEF_sel.m$FPKM<-+2^(BP_genes_DEF_sel.m$value)
              
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_genes_DEF_sel.m_1\n")
                cat(str(BP_genes_DEF_sel.m))
                cat("\n")
                
                
                cat("ensembl_gene_id\n")
                cat(str(BP_genes_DEF_sel.m$ensembl_gene_id))
                cat("\n")
                
                # ########################################
                # quit(status = 1)
              }
              
              
              BP_genes_DEF_sel.m_HET<-droplevels(BP_genes_DEF_sel.m[which(BP_genes_DEF_sel.m$Genotype != "HOM"),])
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("BP_genes_DEF_sel.m_HET_0\n")
                cat(str(BP_genes_DEF_sel.m_HET))
                cat("\n")
                
                
                # quit(status = 1)
                
              }
              
              
              A<-round(summary(BP_genes_DEF_sel.m_HET$value[!is.na(BP_genes_DEF_sel.m_HET$value)]),2)
              
              
              step<-abs(A[6]-A[1])/10
              
              if(Condition_DEBUG == 1)
              {
                cat("summary_BP_genes_DEF_sel.m_HET$value\n")
                cat(sprintf(as.character(names(A))))
                cat("\n")
                cat(sprintf(as.character(A)))
                cat("\n")
              }
              
              if(step == 0)
              {
                
                step<-1
              }
              
              breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
              labels.Rank<-as.character(round(breaks.Rank,3))
              
              levels_ensembl_gene_id<-unique(BP_genes_DEF_sel.m_HET$ensembl_gene_id)
              
              if(Condition_DEBUG == 1)
              {
                cat("labels.Rank:\t")
                cat(sprintf(as.character(labels.Rank)))
                cat("\n")
                
                cat("levels_ensembl_gene_id:\t")
                cat(sprintf(as.character(levels_ensembl_gene_id)))
                cat("\n")
                
                # quit(status = 1)
              }
              
              
              
              
              graph_log2FPKM_BP<-ggplot(data=BP_genes_DEF_sel.m_HET,aes(x=Genotype, y=value, fill=Genotype)) +
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
                scale_y_continuous(name="Gene EXP normalised (log2FPKM)",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
                scale_fill_manual(values=c("gray","#ff1807","#ef8c83"),drop=F)+
                ggeasy::easy_center_title()
              
              # cat("facet_wrap_1:\t")
              # cat("\n")
              
              
              
              graph_log2FPKM_BP <- graph_log2FPKM_BP +
                facet_grid(cols = vars(BP_genes_DEF_sel.m_HET$Cell_Type), scales = "free") +
                theme_cowplot(font_size = 16)+
                theme( strip.background = element_blank(),
                       strip.placement = "inside",
                       strip.text = element_text(size=12, angle=90),
                       panel.spacing = unit(0.2, "lines"), 
                       panel.background=element_rect(fill="white"),
                       panel.border=element_rect(colour="black",size=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
                theme(axis.text.x=element_blank())+
                theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=12))
              
              # if(Condition_DEBUG == 1)
              # {
              #   
              #   cat("graph_log2FPKM_BP_HELLO_WORLD\n")
              #   cat(str(graph_log2FPKM_BP))
              #   cat("\n")
              # }
            }#dim(BP_genes_DEF_sel.m)[1] >0
          }else{
            FLAG_BP<-0}# file.exists("BP_for_LM.rds")
          
          ##### Print per gene ----
          
          path9<-paste(path8,HGNC_sel,'/', sep='')
          
          if (file.exists(path9)){
            
            
            
            
          } else {
            dir.create(file.path(path9))
            
          }
          
          setwd(path9)
          
          if(Condition_DEBUG == 1)
          {
            cat("path9:\t")
            cat(sprintf(as.character(path9)))
            cat("\n")
          }
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_BP:\t")
            cat(sprintf(as.character(FLAG_BP)))
            cat("\n")
            
            # quit(status = 1)
          }
          
        
          graph_DEF<-plot_grid(graph_log2FPKM,graph_adjusted_FPKM,graph_log2FPKM_BP,
                               nrow = 1,
                               ncol = 3,
                               rel_widths=c(0.66,0.66,0.99))
          
        
          
          
          
          
          
          svgname<-paste("DE_EVALUATION_",HGNC_sel,".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= graph_DEF,
                   device="svg",
                   height=10, width=12)
          }
          
          # if(HGNC_sel == "CCR7")
          # {
          #   ######################################################################################
          #   quit(status = 1)
          # }
          
        }#z in 1:length(ALL_TOGETHER_ENSG_array
        # ##############################################################
        # quit(status = 1)
        
      }#length(ALL_TOGETHER_ENSG_array) >0
     
    
    }# i in 1:length(SELECTED_VARS_UPDATED)
    
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
    make_option(c("--SELECTED_VARS_UPDATED"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Tappas_gff"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_INTERVAL_LM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_INTERVAL_LM_Block"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_BP_LM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_BP_LM_Block"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--INTERVAL_GENE_EXP"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_BP"), type="character", default=NULL, 
                metavar="filename", 
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
  
 
  graphical_function(opt)
 
 
  
  
  
}


###########################################################################

system.time( main() )
