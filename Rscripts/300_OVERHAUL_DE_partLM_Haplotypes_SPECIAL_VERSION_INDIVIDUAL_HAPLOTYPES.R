

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

suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)


LM_model = function (option_list)
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
  
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS))
  cat("\n")
  
  #### Haplotypes ----
  
  
  
  Haplotypes = unlist(strsplit(opt$Haplotypes, split=","))
  
  cat("Haplotypes_\n")
  cat(str(Haplotypes))
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
  
  indx.UPDATE<-which(SELECTED_VARS%in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS[-indx.UPDATE]
    
  }else{
    
    SELECTED_VARS_UPDATED = SELECTED_VARS
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
 
  
  ### Read GENES_table----
  
  
  GENES_table<-as.data.frame(fread(file=opt$GENES_table, sep="\t", header=F) , stringsAsFactors=F)
  
  colnames(GENES_table)<-c("chr","start","end","ensembl_gene_id","HGNC")
  
  # cat("GENES_table\n")
  # cat(str(GENES_table))
  # cat("\n")
  
  
  GENES_table_subset<-GENES_table[,c(which(colnames(GENES_table) == "ensembl_gene_id"),
                                     which(colnames(GENES_table) == "HGNC"))]
  
  
  cat("GENES_table_subset\n")
  cat(str(GENES_table_subset))
  cat("\n")
  
  #### Read Thomas matrix of gene expression ----
  
  INTERVAL_GENE_EXP<-as.data.frame(fread(file=opt$INTERVAL_GENE_EXP, sep=",", header=T) , stringsAsFactors=F)
  
  
  
  cat("INTERVAL_GENE_EXP\n")
  cat(str(INTERVAL_GENE_EXP))
  cat("\n")
  # 
  Samples<-unique(INTERVAL_GENE_EXP$sample_id)
  
  # cat("Samples\n")
  # cat(str(Samples))
  # cat("\n")
  
  ENSG_array<-unique(colnames(INTERVAL_GENE_EXP)[-which(colnames(INTERVAL_GENE_EXP) == "sample_id")])
  
  cat("ENSG_array\n")
  cat(str(ENSG_array))
  cat("\n")
  
  # ENSG_array<-ENSG_array[c(1:2)]

  
  #### LOOP ----
  
  Condition_DEBUG <- 0
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    SELECTED_VARS_UPDATED_sel<-SELECTED_VARS_UPDATED
    
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
      
      Condition_DEBUG <- 0
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
        
        
        Haplotype_candidate<-paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__")
        
        FLAG<-Haplotype_candidate[which(Haplotype_candidate%in%Haplotypes)]
        
        if(length(FLAG) >0)
        {
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
              
              array_haplotypes<-levels(Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF$Haplotype)
              
              if(Condition_DEBUG == 1)
              {
                cat("array_haplotypes\n")
                cat(str(array_haplotypes))
                cat("\n")
              }
              
              list_haplotypes<-list()
              Condition_DEBUG <- 0
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
                
                ### Residuals file
                
                setwd(path10)
                
                filename<-paste("DE_LM_HET_RESIDUALS_",comparison_2,'.csv',sep='')
                
                sample_id_vector<-unique(Haplotype_PEER_G_HET_haplotypes_sel$sample_id)
                
                cat("sample_id_vector\n")
                cat(str(sample_id_vector))
                cat("\n")
                
                if (file.exists(filename)){
                  
                  unlink(filename)
                  cat(paste("ensembl_gene_id",paste(sample_id_vector, collapse=","), sep=","),
                      file=filename,sep="\n")
                  
                } else {
                  
                  cat(paste("ensembl_gene_id",paste(sample_id_vector, collapse=","), sep=","),
                      file=filename,sep="\n")
                }
                
                list_RESULT<-list()
                
                
                Condition_DEBUG <- 0
                for(i in 1:length(ENSG_array))
                {
                  ENSG_array_sel<-ENSG_array[i]
                  
                  cat("------------------------------>\t")
                  cat(sprintf(as.character(i)))
                  cat("\t")
                  cat(sprintf(as.character(ENSG_array_sel)))
                  cat("\n")
                  
                  INTERVAL_GENE_EXP_sel<-INTERVAL_GENE_EXP[,c(which(colnames(INTERVAL_GENE_EXP) == "sample_id"),
                                                              which(colnames(INTERVAL_GENE_EXP) %in% ENSG_array_sel))]
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_GENE_EXP_sel_1\n")
                    cat(str(INTERVAL_GENE_EXP_sel))
                    cat("\n")
                  }
                  
                  
                  
                  INTERVAL_GENE_EXP_sel.m<-melt(INTERVAL_GENE_EXP_sel, id.vars=c("sample_id"),value.name = "log2FPKM",variable.name = "ensembl_gene_id")
                  INTERVAL_GENE_EXP_sel.m$FPKM<-+2^(INTERVAL_GENE_EXP_sel.m$log2FPKM)
                  
                  
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
                                                 Haplotype_PEER_G_HET_haplotypes_sel,
                                                 by="sample_id")
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("INTERVAL_GENE_EXP_sel.m_1\n")
                    cat(str(INTERVAL_GENE_EXP_sel.m))
                    cat("\n")
                    cat(sprintf(as.character(colnames(INTERVAL_GENE_EXP_sel.m))))
                    cat("\n")
                    
                    
                    # quit(status = 1)
                    
                  }
                  
                  n_summary<-as.numeric(summary(INTERVAL_GENE_EXP_sel.m$Haplotype))
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("n_summary_0\n")
                    cat(str(n_summary))
                    cat("\n")
                    
                    
                    # quit(status = 1)
                    
                  }
                  
                  n_summary_names<-names(summary(INTERVAL_GENE_EXP_sel.m$Haplotype))
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("n_summary_names_0\n")
                    cat(str(n_summary_names))
                    cat("\n")
                    
                    
                    # quit(status = 1)
                    
                  }
                  
                  RUN_OUT_OF_NAMES<-NULL
                  
                  for(iteration_n_summary_names in 1:length(n_summary_names))
                  {
                    n_summary_names_sel<-n_summary_names[iteration_n_summary_names]
                    n_summary_sel<-n_summary[iteration_n_summary_names]
                    
                    n_summary_string_sel<-paste(n_summary_names_sel,n_summary_sel,sep="__")
                    
                    RUN_OUT_OF_NAMES[iteration_n_summary_names]<-n_summary_string_sel
                    
                  }#iteration_n_summary_names
                  
                  
                  
                  RUN_OUT_OF_NAMES_DEF<-paste(RUN_OUT_OF_NAMES, collapse=";")
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("RUN_OUT_OF_NAMES_DEF_0\n")
                    cat(str(RUN_OUT_OF_NAMES_DEF))
                    cat("\n")
                    
                    
                    # quit(status = 1)
                    
                  }
                  
                  ### FULL linear model of the FPKM (Haplotype+covs) ----
                  
                  unselected_columns<-c("sample_id","ensembl_gene_id","log2FPKM")
                  
                  
                  
                  # Selected_columns<-c("logRatio",
                  #                     "Haplotype",
                  #                     "age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth",
                  #                     "sex",
                  #                     "Season_Winter","Season_Autumn","Season_Spring","Season_Summer",
                  #                     "sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5","sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                  #                     "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                  #                     "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24","PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35",
                  #                     "HCT_PCT___RNA","HGB_g_dL___RNA","IRF_PCT___RNA","MCH_pg___RNA","MCHC_g_dL___RNA","MCV_fL___RNA","RBC_10_12_L___RNA","RDW_SD_fL___RNA","RET_PCT___RNA")
                  
                  INTERVAL_GENE_EXP_sel.m_LM<-INTERVAL_GENE_EXP_sel.m[,-which(colnames(INTERVAL_GENE_EXP_sel.m)%in%unselected_columns)]
                  
                  row.names(INTERVAL_GENE_EXP_sel.m_LM)<-INTERVAL_GENE_EXP_sel.m$sample_id
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("INTERVAL_GENE_EXP_sel.m_LM_0\n")
                    cat(str(INTERVAL_GENE_EXP_sel.m_LM))
                    cat("\n")
                    
                    
                    
                    
                  }
                  
                  
                  model<-lm(FPKM ~ ., data=INTERVAL_GENE_EXP_sel.m_LM)
                  
                  
                  
                  A<-summary(model)
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("A\n")
                    cat(str(A))
                    cat("\n")
                    
                    
                  }
                  
                  
                  results<-A$coefficients
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("results_0\n")
                    cat(str(results))
                    cat("\n")
                    cat(sprintf(as.character(colnames(results))))
                    cat("\n")
                    cat(sprintf(as.character(row.names(results))))
                    cat("\n")
                    cat(sprintf(as.character(results[,4])))
                    cat("\n")
                    
                    results_G<-results[which(row.names(results) == "Haplotype.L"),]
                    
                    cat("results_G\n")
                    cat(str(results_G))
                    cat("\n")
                  }
                  
                  pvalue_Haplotypes_1<-as.numeric(results[which(row.names(results) == "Haplotype.L"),4])
                  coefficient_Haplotypes_1<-as.numeric(results[which(row.names(results) == "Haplotype.L"),1])
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("pvalue_Haplotypes_1\n")
                    cat(str(pvalue_Haplotypes_1))
                    cat("\n")
                    cat("coefficient_Haplotypes_1\n")
                    cat(str(coefficient_Haplotypes_1))
                    cat("\n")
                  }
                  
                  A.df<-as.data.frame(cbind(ENSG_array_sel,pvalue_Haplotypes_1, coefficient_Haplotypes_1,RUN_OUT_OF_NAMES_DEF))
                  
                  colnames(A.df)<-c("ensembl_gene_id","pvalue_Haplotypes_specific_CELL_COUNTS","coefficient_Haplotypes_specific_CELL_COUNTS","n_breakdown_string")
                  
                  A.df$pvalue_Haplotypes_specific_CELL_COUNTS<-as.numeric(A.df$pvalue_Haplotypes_specific_CELL_COUNTS)
                  A.df$coefficient_Haplotypes_specific_CELL_COUNTS<-as.numeric(A.df$coefficient_Haplotypes_specific_CELL_COUNTS)
                  
                  list_RESULT[[i]]<-A.df
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("A.df\n")
                    cat(str(A.df))
                    cat("\n")
                    
                    
                    # quit(status=1)
                  }
                  
                  ### REDUCED linear model of the logRatio (Haplotype+covs) ----
                  
                  unselected_columns<-c("sample_id","ensembl_gene_id","log2FPKM","Haplotype")
                  
                  
                  
                  # Selected_columns<-c("logRatio",
                  #                     "Haplotype",
                  #                     "age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth",
                  #                     "sex",
                  #                     "Season_Winter","Season_Autumn","Season_Spring","Season_Summer",
                  #                     "sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5","sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                  #                     "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                  #                     "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24","PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35",
                  #                     "HCT_PCT___RNA","HGB_g_dL___RNA","IRF_PCT___RNA","MCH_pg___RNA","MCHC_g_dL___RNA","MCV_fL___RNA","RBC_10_12_L___RNA","RDW_SD_fL___RNA","RET_PCT___RNA")
                  
                  INTERVAL_GENE_EXP_sel.m_LM<-INTERVAL_GENE_EXP_sel.m[,-which(colnames(INTERVAL_GENE_EXP_sel.m)%in%unselected_columns)]
                  
                  row.names(INTERVAL_GENE_EXP_sel.m_LM)<-INTERVAL_GENE_EXP_sel.m$sample_id
                  
                  if(Condition_DEBUG == 1)
                  {
                    
                    cat("INTERVAL_GENE_EXP_sel.m_LM_0\n")
                    cat(str(INTERVAL_GENE_EXP_sel.m_LM))
                    cat("\n")
                    
                    
                    
                    
                  }
                  
                  
                  model<-lm(FPKM ~ ., data=INTERVAL_GENE_EXP_sel.m_LM)
                  
                  
                  
                  summary_model<-summary(model)
                  
                  summary_model.m<-melt(summary_model$coefficients)
                  
                  colnames(summary_model.m)[which(colnames(summary_model.m)=="Var1")]<-"Terms"
                  colnames(summary_model.m)[which(colnames(summary_model.m)=="Var2")]<-"Parameters"
                  
                  
                  residual_results=residuals(model)
                  residual_names<-names(residual_results)
                  # colnames(residual_results)<-"residuals_full_model"
                  
                  intercept<-summary_model.m$value[which(summary_model.m$Parameters == "Estimate" & 
                                                           summary_model.m$Terms == "(Intercept)")]
                  
                  
                  
                  # Condition_DEBUG <- 1
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("summary_model.m\n")
                    cat(str(summary_model.m))
                    cat("\n")
                    
                    cat("residual_results_0\n")
                    cat(str(residual_results))
                    cat("\n")
                    
                    cat("residual_names_0\n")
                    cat(str(residual_names))
                    cat("\n")
                    
                    cat("intercept_0\n")
                    cat(str(intercept))
                    cat("\n")
                    
                  }
                  
                  
                  if(length(residual_results) != length(sample_id_vector))
                  {
                    residual_df<-as.data.frame(cbind(residual_names,residual_results),stringsAsFactors=F)
                    colnames(residual_df)<-c("sample_id","residuals")
                    residual_df$residuals<-as.numeric(residual_df$residuals)
                    
                    ### Add intercept to residuals
                    
                    residual_df$residuals<-as.numeric(intercept) + residual_df$residuals
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_ADDED_intercept\n")
                      cat(str(residual_df))
                      cat("\n")
                    }
                    
                    NA_samples<-names(A$na.action)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("NA_samples\n")
                      cat(str(NA_samples))
                      cat("\n")
                    }
                    
                    
                    if(length(NA_samples) >0)
                    {
                      df_NA<- data.frame(matrix(vector(), length(NA_samples), 2,
                                                dimnames=list(c(),
                                                              c("sample_id","residuals"))),
                                         stringsAsFactors=F)
                      
                      df_NA$sample_id<-NA_samples
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("df_NA_0\n")
                        cat(str(df_NA))
                        cat("\n")
                      }
                      
                      df_NA$residuals<-NA
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("df_NA_0.5\n")
                        cat(str(df_NA))
                        cat("\n")
                      }
                      
                      residual_df<-rbind(df_NA,residual_df)
                      
                    }#length(NA_samples) >0
                    
                    
                    
                    residual_df$ensembl_gene_id<-ENSG_array_sel
                    
                    residual_df$sample_id<-factor(residual_df$sample_id,
                                                  levels=sample_id_vector,
                                                  ordered=T)
                    
                    residual_df$residuals<-round(residual_df$residuals,5)
                    
                    residual_df<-residual_df[order(residual_df$sample_id),]
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_1\n")
                      cat(str(residual_df))
                      cat("\n")
                    }
                    
                    residual_df_wide<-as.data.frame(pivot_wider(residual_df,
                                                                id_cols=c("ensembl_gene_id"),
                                                                names_from=sample_id,
                                                                values_from=residuals),
                                                    stringsAsFactors=F)
                    
                    write.table(residual_df_wide,
                                file=filename, append=T,sep=",",
                                quote=F,col.names = F, row.names = F, eol="\n")
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_wide_0\n")
                      cat(str(residual_df_wide))
                      cat("\n")
                      # cat(sprintf(as.character(colnames(residual_df_wide))))
                      # cat("\n")
                      
                      
                    }
                    
                  }else{
                    
                    residual_df<-as.data.frame(cbind(residual_names,residual_results),stringsAsFactors=F)
                    colnames(residual_df)<-c("sample_id","residuals")
                    residual_df$residuals<-as.numeric(residual_df$residuals)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_0\n")
                      cat(str(residual_df))
                      cat("\n")
                    }
                    
                    ### Add intercept to residuals
                    
                    residual_df$residuals<-as.numeric(intercept) + residual_df$residuals
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_ADDED_intercept\n")
                      cat(str(residual_df))
                      cat("\n")
                    }
                    
                    residual_df$ensembl_gene_id<-ENSG_array_sel
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_0.5\n")
                      cat(str(residual_df))
                      cat("\n")
                    }
                    
                    residual_df$sample_id<-factor(residual_df$sample_id,
                                                  levels=sample_id_vector,
                                                  ordered=T)
                    
                    residual_df$residuals<-round(residual_df$residuals,5)
                    
                    
                    residual_df<-residual_df[order(residual_df$sample_id),]
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_1\n")
                      cat(str(residual_df))
                      cat("\n")
                    }
                    
                    residual_df_wide<-as.data.frame(pivot_wider(residual_df,
                                                                id_cols=c("ensembl_gene_id"),
                                                                names_from=sample_id,
                                                                values_from=residuals),
                                                    stringsAsFactors=F)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("residual_df_wide_0\n")
                      cat(str(residual_df_wide))
                      cat("\n")
                    }
                    
                    write.table(residual_df_wide,
                                file=filename, append=T,sep=",",
                                quote=F,col.names = F, row.names = F, eol="\n")
                    
                    
                  }# NA in LM
                }#i in 1:length(ENSG_array)
                
                Condition_DEBUG <- 1
                
                
                if(length(list_RESULT) >0)
                {
                  
                  Results_uncorrected = as.data.frame(data.table::rbindlist(list_RESULT, fill=T), stringsAsFactors=F)
                  
                  Results_uncorrected$comparison<-comparison
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Results_uncorrected_0\n")
                    cat(str(Results_uncorrected))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$ensembl_gene_id)))
                    cat("\n")
                    
                  }
                  
                  Results_uncorrected<-merge(Results_uncorrected,
                                             GENES_table_subset,
                                             by=c("ensembl_gene_id"))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Results_uncorrected_1\n")
                    cat(str(Results_uncorrected))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$ensembl_gene_id)))
                    cat("\n")
                    
                  }
                  
                  list_haplotypes[[l]]<-Results_uncorrected
                  
                  
                  
                }# length(list_RESULT) >0
              }# l in 1:length(array_haplotypes)
              
              if(length(list_haplotypes) >0)
              {
                
                Results_uncorrected_all_comparisons = as.data.frame(data.table::rbindlist(list_haplotypes, fill=T), stringsAsFactors=F)
                
                
                if(Condition_DEBUG == 1)
                {
                  cat("Results_uncorrected_all_comparisons_0\n")
                  cat(str(Results_uncorrected_all_comparisons))
                  cat("\n")
                  cat(str(unique(Results_uncorrected_all_comparisons$ensembl_gene_id)))
                  cat("\n")
                  
                }
                
                path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                
                setwd(path10)
                
                
                saveRDS(Results_uncorrected_all_comparisons,file=paste("DE_LM_Haplotypes_RESULTS_NOMINAL_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'.rds', sep=''))
                
                # ###########################################################
                # quit(status = 1)
                
              }# length(list_haplotypes) >0
              
              
              
            }else{
              
            }#file.exists(filename_rds_haplotype)
            
            
            
          } else {
            
            
          }#file.exists(path10))
        }#length(FLAG) >0
      }#k in 1:length(Proxy_array)
      # ############################################################################################################
      # quit(status = 1)
      
    }#dim(Proxy_file_UPDATED_sel)[1] >0
  }# length(SELECTED_VARS_UPDATED) >0
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
    make_option(c("--GENES_table"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--INTERVAL_GENE_EXP"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SELECTED_VARS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Haplotypes"), type="character", default=NULL,
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
  
  LM_model(opt)
  
  
  
  
}


###########################################################################

system.time( main() )
