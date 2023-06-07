

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

FC_calculation = function(option_list)
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
    
    
    
    ### Read Results_INTERVAL_LM----
    
    
    
    Results_INTERVAL_LM<-as.data.frame(fread(file=opt$Results_INTERVAL_LM, sep="\t", header=T) , stringsAsFactors=T)
    
    
    cat("Results_INTERVAL_LM\n")
    cat(str(Results_INTERVAL_LM))
    cat("\n")
    cat(str(unique(Results_INTERVAL_LM$VAR)))
    cat("\n")
    
  
    
    # quit(status = 1)
    
    
   
    ##### LOOP TO READ ALL VARIABLES -----
    
    
    
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
      
      INTERVAL_covariates_and_PEER_factors_sel<-readRDS(file=paste("INTERVAL_covariates_and_PEER_factors_",SELECTED_VARS_sel,".rds", sep=''))
      
      # cat("INTERVAL_covariates_and_PEER_factors_sel\n")
      # cat(str(INTERVAL_covariates_and_PEER_factors_sel))
      # cat("\n")
      
      residuals_df = as.data.frame(fread(file=paste("DE_LM_HET_RESIDUALS_",SELECTED_VARS_sel,".csv", sep=''),sep=",", header=T), stringsAsFactors=F)
      
      # cat("residuals_df\n")
      # cat(str(residuals_df))
      # cat("\n")
      
      
      
      
     
      
      
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
      
      
   
      ### ONLY FC in INTERVAL
      
      
     
      if(dim(Results_INTERVAL_LM_sel)[1] >0)
      {
        DEF_ENSG_array<-unique(Results_INTERVAL_LM_sel$ensembl_gene_id)
        
        # DEF_ENSG_array<-DEF_ENSG_array[1:10]
        
        if(length(DEF_ENSG_array) >0)
        {
          
          list_FC<-list()
          
          Condition_DEBUG <- 0
          
          for(z in 1:length(DEF_ENSG_array))
          {
            
            DEF_ENSG_array_sel<-DEF_ENSG_array[z]
            
            if(Condition_DEBUG == 1)
            {
              cat("---------------->\t")
              cat(sprintf(as.character(DEF_ENSG_array_sel)))
              cat("\t")
            }
            
            Results_INTERVAL_LM_sel_ENSG_sel<-Results_INTERVAL_LM_sel[which(Results_INTERVAL_LM_sel$ensembl_gene_id%in%DEF_ENSG_array_sel),]
            
            HGNC_sel<-unique(Results_INTERVAL_LM_sel_ENSG_sel$HGNC)
            
            if(Condition_DEBUG == 1)
            {
              cat(sprintf(as.character(HGNC_sel)))
              cat("\n")
            }
            
            if(Condition_DEBUG == 1)
            {
              cat("residuals_df\n")
              cat(str(residuals_df))
              cat("\n")
              
            }
            
            residuals_df_sel_ENSG_sel<-residuals_df[which(residuals_df$ensembl_gene_id == DEF_ENSG_array_sel),]
            
            if(Condition_DEBUG == 1)
            {
              cat("residuals_df_sel_ENSG_sel_\n")
              cat(str(residuals_df_sel_ENSG_sel))
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
            
            
            
            residuals_df_sel_ENSG_sel.m<-melt(residuals_df_sel_ENSG_sel, id.vars=c("ensembl_gene_id"),value.name = "residuals_FPKM",variable.name = "sample_id")
            
            
            if(Condition_DEBUG == 1)
            {
              # cat("residuals_df_sel_ENSG_sel\n")
              # cat(str(residuals_df_sel_ENSG_sel))
              # cat("\n")
              
              cat("residuals_df_sel_ENSG_sel.m\n")
              cat(str(residuals_df_sel_ENSG_sel.m))
              cat("\n")
            }
            
           
            #### Merge with covariates and Keep to HET ----
            
            
            
            
            residuals_df_sel_ENSG_sel.m<-merge(residuals_df_sel_ENSG_sel.m,
                                           INTERVAL_covariates_and_PEER_factors_sel_subset,
                                           by=c("sample_id"))
            A<-round(summary(residuals_df_sel_ENSG_sel.m$residuals_FPKM[!is.na(residuals_df_sel_ENSG_sel.m$residuals_FPKM)]),2)
            # A2<-round(summary(residuals_df_sel_ENSG_sel.m$FPKM[!is.na(residuals_df_sel_ENSG_sel.m$FPKM)]),2)
            
            residuals_df_sel_ENSG_sel.m$residuals_FPKM_no_negative<-residuals_df_sel_ENSG_sel.m$residuals_FPKM+abs(A[1])
            
            if(Condition_DEBUG == 1)
            {
              
              cat("residuals_df_sel_ENSG_sel.m_3\n")
              cat(str(residuals_df_sel_ENSG_sel.m))
              cat("\n")
              
              # quit(status = 1)
              
            }
            
            residuals_df_sel_ENSG_sel.m_HET<-droplevels(residuals_df_sel_ENSG_sel.m[(as.numeric(residuals_df_sel_ENSG_sel.m$Genotype) <= 2),])
            
            if(Condition_DEBUG == 1)
            {
              
              cat("residuals_df_sel_ENSG_sel.m_HET_0\n")
              cat(str(residuals_df_sel_ENSG_sel.m_HET))
              cat("\n")
              
              
              # quit(status = 1)
              
            }
            
           
            
            residuals_df_sel_ENSG_sel.m_HET.dt<-data.table(residuals_df_sel_ENSG_sel.m_HET, key=c("ensembl_gene_id","Genotype"))
            
            Summary_table<-as.data.frame(residuals_df_sel_ENSG_sel.m_HET.dt[,.(Median_residuals_FPKM_no_negative =median(residuals_FPKM_no_negative)
            ),by=key(residuals_df_sel_ENSG_sel.m_HET.dt)], stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_LogFC_0\n")
              cat(str(Summary_table))
              cat("\n")
              
            }
            
            Summary_table_wide<-as.data.frame(pivot_wider(Summary_table,
                                                          id_cols=ensembl_gene_id,
                                                          names_from=Genotype,
                                                          values_from=Median_residuals_FPKM_no_negative), stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_wide_LogFC_0\n")
              cat(str(Summary_table_wide))
              cat("\n")
              
            }
            
            Results_INTERVAL_LM_sel_ENSG_sel$Median_residuals_FPKM_no_negative_HOM_REF<-Summary_table_wide$HOM_REF
            Results_INTERVAL_LM_sel_ENSG_sel$Median_residuals_FPKM_no_negative_HET<-Summary_table_wide$HET
            Results_INTERVAL_LM_sel_ENSG_sel$LogFC<-foldchange2logratio(foldchange(Summary_table_wide$HET,Summary_table_wide$HOM_REF))
            
            
            if(Condition_DEBUG == 1)
            {
              cat("Results_INTERVAL_LM_sel_ENSG_sel_for_dot_plot_FOR_LOGFC_1:\t")
              cat(str(Results_INTERVAL_LM_sel_ENSG_sel))
              cat("\n")
              
            }
            
            indx.int<-c(which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "VAR"),
                        which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "ensembl_gene_id"),which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "HGNC"),
                        which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "Genome_wide_minuslogpvalue"),which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "Genome_wide_Beta"),
                        which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "Median_residuals_FPKM_no_negative_HOM_REF"),which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "Median_residuals_FPKM_no_negative_HET"),
                        which(colnames(Results_INTERVAL_LM_sel_ENSG_sel) == "LogFC"))

            Results_INTERVAL_LM_sel_ENSG_sel_subset<-unique(Results_INTERVAL_LM_sel_ENSG_sel[,indx.int])
            
            if(Condition_DEBUG == 1)
            {
              cat("Results_INTERVAL_LM_sel_ENSG_sel_subset:\t")
              cat(str(Results_INTERVAL_LM_sel_ENSG_sel_subset))
              cat("\n")
              
            }
            
            list_FC[[z]]<-Results_INTERVAL_LM_sel_ENSG_sel_subset
            
            # ##############################
            # quit(status = 1)
          
          }#z in 1:length(DEF_ENSG_array)
          
          Condition_DEBUG <- 0
          
          if(length(list_FC) >0)
          {
            
            FC_per_VAR = as.data.frame(data.table::rbindlist(list_FC, fill=T), stringsAsFactors=F)
            
            
            if(Condition_DEBUG == 1)
            {
              cat("FC_per_VAR_0\n")
              cat(str(FC_per_VAR))
              cat("\n")
              cat(str(unique(FC_per_VAR$VAR)))
              cat("\n")
              cat(str(unique(FC_per_VAR$ensembl_gene_id)))
              cat("\n")
            }
            
            setwd(path6)
            
            saveRDS(file="DE_RESULTS_FOR_GSEA.rds",FC_per_VAR)
            
            # #############################
            # quit(status = 1)
            
          }# length(list_FC) >0
          
        }#length(DEF_ENSG_array) >0
        
      }# INTERVAL exist
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
    make_option(c("--Results_INTERVAL_LM"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_INTERVAL_LM_Block"), type="character", default=NULL,
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
  
  FC_calculation(opt)
  
 
 
  
  
  
}


###########################################################################

system.time( main() )
