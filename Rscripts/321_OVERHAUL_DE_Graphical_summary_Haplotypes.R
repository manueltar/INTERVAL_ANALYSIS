

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
library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")


opt = NULL

# options(warn = 1)


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
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS_UPDATED  = unlist(strsplit(opt$SELECTED_VARS_UPDATED , split=","))
  
  # SELECTED_VARS_UPDATED<-"chr5_35476470_G_T"
  
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
  
  
  ##### LOOP TO READ ALL VARIABLES -----
  
  
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
    
    
    
    
    ### Read initial selection----
    
    
    Supp4_CURATED_WITH_PHENOTYPES<-as.data.frame(readRDS(file=opt$Supp4_CURATED_WITH_PHENOTYPES) , stringsAsFactors=F)
    
    cat("Supp4_CURATED_WITH_PHENOTYPES_0\n")
    cat(str(Supp4_CURATED_WITH_PHENOTYPES))
    cat("\n")
    cat(str(unique(Supp4_CURATED_WITH_PHENOTYPES$VAR)))
    cat("\n")
    
    
    
    # quit(status=1)
    
    
    ### Read Results_INTERVAL_LM----
    
    
    
    Results_INTERVAL_LM<-as.data.frame(fread(file=opt$Results_INTERVAL_LM, sep="\t", header=T) , stringsAsFactors=T)
    
    
    cat("Results_INTERVAL_LM\n")
    cat(str(Results_INTERVAL_LM))
    cat("\n")
    cat(str(unique(Results_INTERVAL_LM$VAR)))
    cat("\n")
    
   
    # quit(status = 1)
    
    
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
    
   
    
    breaks.size<-sort(c(1,seq(0,8, by=2)))
    labels.size<-as.character(breaks.size)
    
    cat("labels.size\n")
    cat(str(labels.size))
    cat("\n")
    
    # SELECTED_VARS_UPDATED<-"chr16_67250992_C_T"
    
    for(i in 1:length(SELECTED_VARS_UPDATED))
    {
      SELECTED_VARS_UPDATED_sel<-SELECTED_VARS_UPDATED[i]
      
      # cat("MAIN_VAR---------->\t")
      # cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
      # cat("\n")
      
      
      cat("MAIN_VAR---------->\t")
      cat(sprintf(as.character(SELECTED_VARS_UPDATED_sel)))
      cat("\n")
      
      Proxy_file_UPDATED_sel<-Proxy_file_UPDATED[which(Proxy_file_UPDATED$VAR%in%SELECTED_VARS_UPDATED_sel),]
      
      
      if(dim(Proxy_file_UPDATED_sel)[1] >0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("Proxy_file_UPDATED_sel\n")
          cat(str(Proxy_file_UPDATED_sel))
          cat("\n")
          cat(str(unique(Proxy_file_UPDATED_sel$VAR)))
          cat("\n")
          cat(str(unique(Proxy_file_UPDATED_sel$Proxy_VAR)))
          cat("\n")
        }
        
        Proxy_array<-unique(Proxy_file_UPDATED_sel$Proxy_VAR)
        
        if(Condition_DEBUG == 1)
        {
          cat("Proxy_array\n")
          cat(str(Proxy_array))
          cat("\n")
        }
        
        # Proxy_array<-"chr16_67520435_A_G"
        
        # Proxy_array<-"chr5_35396231_C_T"
        
        for(k in 1:length(Proxy_array))
        {
          Proxy_array_sel<-Proxy_array[k]
          
          
          
          cat("Proxy_VAR---------->\t")
          cat(sprintf(as.character(Proxy_array_sel)))
          cat("\n")
          
          
          Proxy_df<-Proxy_file_UPDATED_sel[which(Proxy_file_UPDATED_sel$Proxy_VAR == Proxy_array_sel),]
          
          if(Condition_DEBUG == 1)
          {
            cat("Proxy_df\n")
            cat(str(Proxy_df))
            cat("\n")
            cat(str(unique(Proxy_df$VAR)))
            cat("\n")
            cat(str(unique(Proxy_df$Proxy_VAR)))
            cat("\n")
          }
          
          path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
          
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path10)){
            
            setwd(path10)
            
            filename_rds_haplotype<-paste("INTERVAL_covariates_and_PEER_factors_Haplotype_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),".rds", sep='')
            
            if (file.exists(filename_rds_haplotype)){
              
              Haplotype_PEER_G<-readRDS(file=filename_rds_haplotype)
              
              if(Condition_DEBUG == 1)
              {
                cat("Haplotype_PEER_G\n")
                cat(str(Haplotype_PEER_G))
                cat("\n")
              }
              
              Haplotype_PEER_G_HET_haplotypes<-droplevels(Haplotype_PEER_G[(as.numeric(Haplotype_PEER_G$Haplotype) <= 4),])
              
              Haplotype_PEER_G_HET_haplotypes$Haplotype<-factor(Haplotype_PEER_G_HET_haplotypes$Haplotype,
                                                                levels=c('HOM_REF','HET|HOM_REF','HOM_REF|HET','HET|HET'),
                                                                ordered=T)
              
              if(Condition_DEBUG == 1)
              {
                cat("Haplotype_PEER_G_HET_haplotypes\n")
                cat(str(Haplotype_PEER_G_HET_haplotypes))
                cat("\n")
              }
              
              
              
              Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF<-droplevels(Haplotype_PEER_G_HET_haplotypes[(as.numeric(Haplotype_PEER_G_HET_haplotypes$Haplotype) > 1),])
              
              
              if(Condition_DEBUG == 1)
              {
                cat("Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF\n")
                cat(str(Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF))
                cat("\n")
              }
              
              Haplotype_PEER_G_HET_haplotypes_HOM_REF<-droplevels(Haplotype_PEER_G_HET_haplotypes[(as.numeric(Haplotype_PEER_G_HET_haplotypes$Haplotype) == 1),])
              
              if(Condition_DEBUG == 1)
              {
                cat("Haplotype_PEER_G_HET_haplotypes_HOM_REF\n")
                cat(str(Haplotype_PEER_G_HET_haplotypes_HOM_REF))
                cat("\n")
              }
              
              samples_HOM_REF<-Haplotype_PEER_G_HET_haplotypes_HOM_REF$sample_id
              
              if(Condition_DEBUG == 1)
              {
                cat("samples_HOM_REF\n")
                cat(str(samples_HOM_REF))
                cat("\n")
              }
              
              Results_INTERVAL_LM_sel<-Results_INTERVAL_LM[which(Results_INTERVAL_LM$VAR == SELECTED_VARS_UPDATED_sel &
                                                                   Results_INTERVAL_LM$Proxy_VAR == Proxy_array_sel),]
              
              if(Condition_DEBUG == 1)
              {
                cat("Results_INTERVAL_LM_sel\n")
                cat(str(Results_INTERVAL_LM_sel))
                cat("\n")
                cat(str(unique(Results_INTERVAL_LM_sel$VAR)))
                cat("\n")
                cat(str(unique(Results_INTERVAL_LM_sel$Proxy_VAR)))
                cat("\n")
              }
              
              ENSG_array<-unique(c(Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$CIS_gene_minuslogpvalue)],
                                   Results_INTERVAL_LM_sel$ensembl_gene_id[!is.na(Results_INTERVAL_LM_sel$Block_PCHiC_minuslogpvalue)]))
              
              # ENSG_array<-"RP11-297D21.4"
              
              if(length(ENSG_array) >0)
              {
                if(Condition_DEBUG == 1)
                {
                  cat("ENSG_array_\n")
                  cat(str(ENSG_array))
                  cat("\n")
                }
                
                array_haplotypes<-levels(Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF$Haplotype)
                
                if(Condition_DEBUG == 1)
                {
                  cat("array_haplotypes\n")
                  cat(str(array_haplotypes))
                  cat("\n")
                }
                
                list_haplotypes_residuals<-list()
                
                for(l in 1:length(array_haplotypes))
                {
                  
                  array_haplotypes_sel<-array_haplotypes[l]
                  
                  
                  cat("---------->\t")
                  cat(sprintf(as.character(array_haplotypes_sel)))
                  cat("\n")
                  
                  Haplotype_PEER_G_HET_haplotypes_sel<-droplevels(Haplotype_PEER_G_HET_haplotypes[which(Haplotype_PEER_G_HET_haplotypes$Haplotype%in%c("HOM_REF",array_haplotypes_sel)),])
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Haplotype_PEER_G_HET_haplotypes_sel\n")
                    cat(str(Haplotype_PEER_G_HET_haplotypes_sel))
                    cat("\n")
                  }
                  
                  comparison<-paste(levels(Haplotype_PEER_G_HET_haplotypes_sel$Haplotype), collapse = "__")
                  
                  comparison_2<-gsub("\\|","_",comparison)
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("comparison\n")
                    cat(str(comparison))
                    cat("\n")
                    
                    cat("comparison_2\n")
                    cat(str(comparison_2))
                    cat("\n")
                  }
                  
                  
                  
                  if (file.exists(paste("DE_LM_HET_RESIDUALS_",comparison_2,'.csv',sep='')))
                  {
                    
                    residuals_df = read.csv(file=paste("DE_LM_HET_RESIDUALS_",comparison_2,'.csv',sep=''))
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residuals_df\n")
                      cat(str(residuals_df))
                      cat("\n")
                    }
                    
                    residuals_df_sel<-residuals_df[which(residuals_df$ensembl_gene_id%in%ENSG_array),]
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residuals_df_sel\n")
                      cat(str(residuals_df_sel))
                      cat("\n")
                    }
                    
                    if(dim(residuals_df_sel)[1] >0)
                    {
                      
                      if(l ==1)
                      {
                        residuals_df_sel.m<-melt(residuals_df_sel, id.vars=c("ensembl_gene_id"),value.name = "residuals_FPKM",variable.name = "sample_id")
                        
                        
                      }#l >1
                      else{
                        
                        residuals_df_sel<-residuals_df_sel[,-which(colnames(residuals_df_sel)%in%samples_HOM_REF)]
                        
                        residuals_df_sel.m<-melt(residuals_df_sel, id.vars=c("ensembl_gene_id"),value.name = "residuals_FPKM",variable.name = "sample_id")
                        
                        
                      }
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("residuals_df_sel.m\n")
                        cat(str(residuals_df_sel.m))
                        cat("\n")
                      }
                      
                      list_haplotypes_residuals[[l]]<-residuals_df_sel.m
                    }# dim(residuals_df_sel)[1] >0 check chr5_35476470_G_T/Haplotypes/chr5_35476470_G_T__chr5_35396231_C_T/
                  }#file.exists(paste("DTU_LM_HET_RESIDUALS_",comparison_2,'.csv',sep='')))
                }#l in 1:length(array_haplotypes)
                
                if(length(list_haplotypes_residuals) >0)
                {
                  
                  residuals_df_ALL_Haplotypes = unique(as.data.frame(data.table::rbindlist(list_haplotypes_residuals, fill=T), stringsAsFactors=F))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("residuals_df_ALL_Haplotypes_0\n")
                    cat(str(residuals_df_ALL_Haplotypes))
                    cat("\n")
                    
                  
                    
                    #quit(status = 1)
                  }
                  
                }# length(list_haplotypes_residuals) >0
                
                
               
                Results_INTERVAL_LM_sel_ENSG_sel<-Results_INTERVAL_LM_sel[which(Results_INTERVAL_LM_sel$ensembl_gene_id%in%ENSG_array),]
                
                if(dim(Results_INTERVAL_LM_sel_ENSG_sel)[1] >0)
                {
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
                  
                  DE_gene_df<-band.genes<-data.frame(matrix(ncol=6,nrow=dim(Results_INTERVAL_LM_sel_ENSG_sel)[1], 
                                                            dimnames=list(NULL, c("VAR", "Proxy_VAR","comparison",
                                                                                  "ensembl_gene_id",
                                                                                  "Significance",
                                                                                  "RNASeq_source"))),
                                                     stringsAsFactors = F)
                  
                  DE_gene_df$VAR<-Results_INTERVAL_LM_sel_ENSG_sel$VAR
                  DE_gene_df$Proxy_VAR<-Results_INTERVAL_LM_sel_ENSG_sel$Proxy_VAR
                  DE_gene_df$comparison<-Results_INTERVAL_LM_sel_ENSG_sel$comparison
                  DE_gene_df$ensembl_gene_id<-Results_INTERVAL_LM_sel_ENSG_sel$ensembl_gene_id
                  
                  
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
                    Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "comparison"),
                                                                         which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                         which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                         which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                                         which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"),
                                                                         which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta_Z_score"))]

                    colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Haplotypes"
                    colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Haplotypes"
                    colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Haplotypes_Z_score"
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Block_genes_INTERVAL_subset_0\n")
                      cat(str(Block_genes_INTERVAL_subset))
                      cat("\n")
                      
                    }
                    
                    CIS_gene_INTERVAL_subset<-CIS_gene_INTERVAL[,c(which(colnames(CIS_gene_INTERVAL) == "comparison"),
                                                                   which(colnames(CIS_gene_INTERVAL) == "ensembl_gene_id"),
                                                                   which(colnames(CIS_gene_INTERVAL) == "HGNC"),
                                                                   which(colnames(CIS_gene_INTERVAL) == "CIS_gene_minuslogpvalue"),
                                                                   which(colnames(CIS_gene_INTERVAL) == "CIS_Beta"),
                                                                   which(colnames(CIS_gene_INTERVAL) == "CIS_Beta_Z_score"))]
                    
                    colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_gene_minuslogpvalue")]<-"ajusted.minuslogpvalue_Haplotypes"
                    colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta")]<-"coefficient_Haplotypes"
                    colnames(CIS_gene_INTERVAL_subset)[which(colnames(CIS_gene_INTERVAL_subset)== "CIS_Beta_Z_score")]<-"coefficient_Haplotypes_Z_score"
                    
                    
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
                                      by=c("comparison","ensembl_gene_id"),
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
                    
                    Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "comparison"),
                                                                         which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                         which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                         which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_minuslogpvalue"),
                                                                         which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta"),
                                                                         which(colnames(Block_genes_INTERVAL) == "Block_PCHiC_Beta_Z_score"))]
                    
                    colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_minuslogpvalue")]<-"ajusted.minuslogpvalue_Haplotypes"
                    colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta")]<-"coefficient_Haplotypes"
                    colnames(Block_genes_INTERVAL_subset)[which(colnames(Block_genes_INTERVAL_subset)== "Block_PCHiC_Beta_Z_score")]<-"coefficient_Haplotypes_Z_score"
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Block_genes_INTERVAL_subset_0\n")
                      cat(str(Block_genes_INTERVAL_subset))
                      cat("\n")
                      
                    }
                    
                    DE_gene_df<-merge(DE_gene_df,
                                      Block_genes_INTERVAL_subset,
                                      by=c("comparison","ensembl_gene_id"),
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
                  
                  DE_gene_df$Significance[which(DE_gene_df$ajusted.minuslogpvalue_Haplotypes >= 1.3)]<-"YES"
                  DE_gene_df$Significance[which(DE_gene_df$ajusted.minuslogpvalue_Haplotypes < 1.3)]<-"NO"
                  
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
                  
                  
                  ### update break.size
                  
                  local_minuslogpval_max<-max(DE_gene_df$ajusted.minuslogpvalue_Haplotypes[!is.na(DE_gene_df$ajusted.minuslogpvalue_Haplotypes)])
                  
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
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("breaks.size_updated\n")
                      cat(sprintf(as.character(breaks.size_updated)))
                      cat("\n")
                    }
                    
                    # quit(status = 1)
                  }
                  
                  labels.size_updated<-as.character(breaks.size_updated)
                  
                  
                  
                  
                  
                  ### Use Beta from the full model of logRatio ---
                  
                  ### Labels beta
                  
                  indx.finite<-is.finite(DE_gene_df$coefficient_Haplotypes_Z_score)
                  
                  check.finite<-sum(indx.finite)
                  
                  if(check.finite >0)
                  {
                    A<-summary(DE_gene_df$coefficient_Haplotypes_Z_score[indx.finite])
                    
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
                    
                    breaks.Beta_INTERVAL<-sort(unique(round(c(0,candidate_vector,-0.25,0.25),3)))
                    labels.Beta_INTERVAL<-as.character(breaks.Beta_INTERVAL)
                    
                  }else{
                    
                    breaks.Beta_INTERVAL<-sort(unique(round(c(0,-0.25,0.25),3)))
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
                  
                  
                  ### Graph ----
                  
                  DE_gene_df$comparison<-factor(DE_gene_df$comparison,
                                                levels=c('HOM_REF__HET|HOM_REF',
                                                         'HOM_REF__HOM_REF|HET',
                                                         'HOM_REF__HET|HET'),
                                                ordered=T)
                  
                  DE_gene_df<-DE_gene_df[order(DE_gene_df$ensembl_gene_id),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("DE_gene_df_for_dot_plot_2:\t")
                    cat(str(DE_gene_df))
                    cat("\n")
                    cat(sprintf(as.character(names(summary(DE_gene_df$comparison)))))
                    cat("\n")
                    cat(sprintf(as.character(summary(DE_gene_df$comparison))))
                    cat("\n")
                  }
                  
                  dotplot_INTERVAL<-ggplot(data=DE_gene_df,
                                           aes(y=HGNC,
                                               x=RNASeq_source)) +
                    geom_point(aes(color=Significance,
                                   fill=coefficient_Haplotypes_Z_score,
                                   size=ajusted.minuslogpvalue_Haplotypes), stroke=1, shape=21)+
                    scale_color_manual(values=c("gray","black"),name='Significant', drop=F)+
                    scale_size(range = c(0,20), name='-log10pval',
                               breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
                    scale_fill_gradient2(
                      low = "blue", 
                      mid = "white", 
                      high = "red", 
                      midpoint = 0,
                      breaks=breaks.Beta_INTERVAL,labels=labels.Beta_INTERVAL,
                      limits=c(breaks.Beta_INTERVAL[1]-0.01,breaks.Beta_INTERVAL[length(breaks.Beta_INTERVAL)]+0.01),name=paste('Effect size','Z-score', sep="\n"),na.value = "gray")+
                    scale_y_discrete(name=NULL, drop=F)+
                    theme_classic()+
                    ggeasy::easy_center_title()
                  
                  dotplot_INTERVAL<-dotplot_INTERVAL+
                    facet_grid(cols = vars(DE_gene_df$comparison), scales='free_x', space='free_x') +
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
                  
                  
                  
                  graph_DEF<-plot_grid(dotplot_INTERVAL,
                                       nrow = 1,
                                       ncol = 1)
                  
                  path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
                  path8<-paste(path7,SELECTED_VARS_UPDATED_sel,'/', sep='')
                  path9<-paste(path8,'Haplotypes','/', sep='')
                  if (file.exists(path9)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path9))
                    
                  }
                  
                  path10<-paste(path9,paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                  
                  if (file.exists(path10)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path10))
                    
                  }
                  
                  ##### PRINT DOT PLOTS -------------------------
                  
                  setwd(path10)
                  
                  svgname<-paste("DE_Series_1_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),".svg",sep='')
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
                  
                  ##### PRINT Per GENE violin plots -------------------------
                  
                  for(z in 1:length(ENSG_array))
                  {
                    
                    ENSG_array_sel<-ENSG_array[z]
                    
                    cat("---------------->\t")
                    cat(sprintf(as.character(ENSG_array_sel)))
                    cat("\t")
                    
                    DE_gene_df_ENSG_sel<-DE_gene_df[which(DE_gene_df$ensembl_gene_id%in%ENSG_array_sel),]
                    
                    HGNC_sel<-unique(DE_gene_df_ENSG_sel$HGNC)
                    
                    cat(sprintf(as.character(HGNC_sel)))
                    cat("\n")
                    
                    
                    Haplotype_PEER_G_HET_haplotypes_subset<-Haplotype_PEER_G_HET_haplotypes[,c(which(colnames(Haplotype_PEER_G_HET_haplotypes) == "Haplotype"),
                                                                                               which(colnames(Haplotype_PEER_G_HET_haplotypes) == "sample_id"))]
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Haplotype_PEER_G_HET_haplotypes_subset_\n")
                      cat(str(Haplotype_PEER_G_HET_haplotypes_subset))
                      cat("\n")
                      
                    }
                    
                    INTERVAL_GENE_EXP_sel<-INTERVAL_GENE_EXP[,c(which(colnames(INTERVAL_GENE_EXP) == "sample_id"),
                                                                which(colnames(INTERVAL_GENE_EXP) %in% ENSG_array_sel))]
                    if(Condition_DEBUG == 1)
                    {
                      cat("INTERVAL_GENE_EXP_sel_1\n")
                      cat(str(INTERVAL_GENE_EXP_sel))
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
                    
                    if(length(list_haplotypes_residuals) >0)
                    {
                      INTERVAL_GENE_EXP_sel.m<-merge(INTERVAL_GENE_EXP_sel.m,
                                                     residuals_df_ALL_Haplotypes,
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
                                                     Haplotype_PEER_G_HET_haplotypes_subset,
                                                     by=c("sample_id"))
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("INTERVAL_GENE_EXP_sel.m_3\n")
                        cat(str(INTERVAL_GENE_EXP_sel.m))
                        cat("\n")
                        
                        # quit(status = 1)
                        
                      }
                      
                      INTERVAL_GENE_EXP_sel.m_REDUCED<-INTERVAL_GENE_EXP_sel.m[(as.numeric(INTERVAL_GENE_EXP_sel.m$Haplotype) <= 4),]
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("INTERVAL_GENE_EXP_sel.m_REDUCED_0\n")
                        cat(str(INTERVAL_GENE_EXP_sel.m_REDUCED))
                        cat("\n")
                        
                        
                        # quit(status = 1)
                        
                      }
                      
                      ##### Violin visual log2FPKM unadjusted----
                      
                      A<-round(summary(INTERVAL_GENE_EXP_sel.m_REDUCED$log2FPKM[!is.na(INTERVAL_GENE_EXP_sel.m_REDUCED$log2FPKM)]),2)
                      
                      
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
                      
                      levels_ensembl_gene_id<-unique(INTERVAL_GENE_EXP_sel.m_REDUCED$ensembl_gene_id)
                      
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
                      
                      
                      
                      
                      graph_log2FPKM<-ggplot(data=INTERVAL_GENE_EXP_sel.m_REDUCED,aes(x=Haplotype, y=log2FPKM, fill=Haplotype)) +
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
                        scale_fill_brewer(palette = "Dark2", drop=F)+
                        ggtitle(paste(unique(INTERVAL_GENE_EXP_sel.m_REDUCED$ensembl_gene_id),unique(INTERVAL_GENE_EXP_sel.m_REDUCED$HGNC),sep=' '))+
                        ggeasy::easy_center_title()
                      
                      
                      ##### Violin visual FPKM_adjusted reduced model ----
                      
                      INTERVAL_GENE_EXP_sel.m_REDUCED$Log2_residuals_FPKM_no_negative<-log2(INTERVAL_GENE_EXP_sel.m_REDUCED$residuals_FPKM_no_negative)
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("INTERVAL_GENE_EXP_sel.m_REDUCED_FPKM_reduced_model\n")
                        cat(str(INTERVAL_GENE_EXP_sel.m_REDUCED))
                        cat("\n")
                        
                        
                        # quit(status = 1)
                        
                      }
                      
                      indx.infinite<-is.infinite(INTERVAL_GENE_EXP_sel.m_REDUCED$Log2_residuals_FPKM_no_negative[!is.na(INTERVAL_GENE_EXP_sel.m_REDUCED$Log2_residuals_FPKM_no_negative)])
                      
                      check.infinite<-sum(indx.infinite)
                      
                      if(check.infinite ==0)
                      {
                        A<-round(summary(INTERVAL_GENE_EXP_sel.m_REDUCED$Log2_residuals_FPKM_no_negative[!is.na(INTERVAL_GENE_EXP_sel.m_REDUCED$Log2_residuals_FPKM_no_negative)]),3)
                        
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("INTERVAL_GENE_EXP_sel.m_REDUCED$Log2_residuals_FPKM_no_negative\n")
                          cat(sprintf(as.character(names(A))))
                          cat("\n")
                          cat(sprintf(as.character(A)))
                          cat("\n")
                        }
                        
                        step<-abs(A[6]-A[1])/10
                        
                        if(step == 0)
                        {
                          
                          step<-1
                        }
                        
                        breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
                        labels.Rank<-as.character(round(breaks.Rank,3))
                        
                      }else{
                        
                        breaks.Rank<-unique(seq(from= -4, to=4,by=0.5))
                        labels.Rank<-as.character(round(breaks.Rank,3))
                        
                        
                      }# length(indx.finite) >0
                      
                      
                      
                      
                      levels_ensembl_gene_id<-unique(INTERVAL_GENE_EXP_sel.m_REDUCED$ensembl_gene_id)
                      
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
                      
                      
                      
                      
                      graph_adjusted_FPKM<-ggplot(data=INTERVAL_GENE_EXP_sel.m_REDUCED,aes(x=Haplotype, y=Log2_residuals_FPKM_no_negative, fill=Haplotype)) +
                        geom_violin()+
                        stat_summary(fun = median, fun.min = median, fun.max = median,
                                     geom = "crossbar", width = 0.5)+
                        theme_bw()+
                        theme(axis.title.y=element_text(size=24, family="sans"),
                              axis.title.x=element_text(size=24, family="sans"),
                              axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                              axis.text.x=element_text(angle=45,size=18,hjust=1,vjust=1, color="black", family="sans"),
                              legend.title=element_text(size=16,color="black", family="sans"),
                              legend.text=element_text(size=12,color="black", family="sans"))+
                        scale_x_discrete(name=NULL, drop=F)+
                        scale_y_continuous(name="Residuals model Gene EXP log2FPKM adjusted for covariates",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
                        scale_fill_brewer(palette = "Dark2", drop=F)+
                        ggtitle(paste(unique(INTERVAL_GENE_EXP_sel.m_REDUCED$ensembl_gene_id),unique(INTERVAL_GENE_EXP_sel.m_REDUCED$HGNC),sep=' '))+
                        ggeasy::easy_center_title()
                      
                      ##### Print per gene ----
                      path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
                      
                      if (file.exists(path7)){
                        
                        
                        
                        
                      } else {
                        dir.create(file.path(path7))
                        
                      }
                      
                      path8<-paste(path7,SELECTED_VARS_UPDATED_sel,'/', sep='')
                      
                      if (file.exists(path8)){
                        
                        
                        
                        
                      } else {
                        dir.create(file.path(path8))
                        
                      }
                      
                      path9<-paste(path8,'Haplotypes','/', sep='')
                      
                      if (file.exists(path9)){
                        
                        
                        
                        
                      } else {
                        dir.create(file.path(path9))
                        
                      }
                      
                      path10<-paste(path9,paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                      
                      if (file.exists(path10)){
                        
                        
                        
                        
                      } else {
                        dir.create(file.path(path10))
                        
                      }
                      
                      path11<-paste(path10,HGNC_sel,'/', sep='')
                      
                      if (file.exists(path11)){
                        
                        
                        
                        
                      } else {
                        dir.create(file.path(path11))
                        
                      }
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("path11:\t")
                        cat(sprintf(as.character(path11)))
                        cat("\n")
                      }
                      
                      setwd(path11)
                      
                      
                      
                      
                      graph_DEF<-plot_grid(graph_adjusted_FPKM,
                                           nrow = 1,
                                           ncol = 1)
                      
                      
                      
                      
                      svgname<-paste("DE_EVALUATION_",HGNC_sel,".svg",sep='')
                      makesvg = TRUE
                      
                      if (makesvg == TRUE)
                      {
                        ggsave(svgname, plot= graph_DEF,
                               device="svg",
                               height=10, width=12)
                      }
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("svgname:\t")
                        cat(sprintf(as.character(svgname)))
                        cat("\n")
                      }
                      
                    }else{
                      
                      cat("WARNING!! ABSENT GENE FROMN RESIDUALS\n")
                    }# Possible error of hard drive capacity length(list_haplotypes_residuals) >0
                   
                    
                    # ###################################################################
                    # quit(status = 1)
                    
                  }#z in 1:length(ENSG_array

                }#dim(Results_INTERVAL_LM_sel_ENSG_sel)[1] >0
              }#length(ENSG_array) >0
            }#file.exists(filename_rds_haplotype)
          }#file.exists(path10)
        }#k in 1:length(Proxy_array)
      }#dim(Proxy_file_UPDATED_sel)[1] >0
    }#i in 1:length(SELECTED_VARS_UPDATED)
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
    make_option(c("--Supp4_CURATED_WITH_PHENOTYPES"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--INTERVAL_GENE_EXP"), type="character", default=NULL,
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
  
 
  graphical_function(opt)
 
 
  
  
  
}


###########################################################################

system.time( main() )
