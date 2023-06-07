

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
  
 
  ##### LOOP TO READ ALL VARIABLES -----
  
  
  Condition_DEBUG <- 0
  
  
  
  if(length(SELECTED_VARS_UPDATED) >0)
  {
    list_GW_DEF<-list()
    list_Block_DEF<-list()
    list_CIS_DEF<-list()
    
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
        
        list_GW<-list()
        list_Block<-list()
        list_CIS<-list()
        
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
            ### DTU genome wide
            
            filename="DTU_genome_wide.rds"
            
            n_SIG_GW<-0
            n_tested_GW<-0
            Perc_SIG_GW<-NA
            
            if (file.exists(filename)) {
              DTU_genome_wide = readRDS(filename)
              
              if(Condition_DEBUG == 1)
              {
                cat("DTU_genome_wide_0\n")
                cat(str(DTU_genome_wide))
                cat("\n")
              }
              
              DTU_genome_wide_SIG<-DTU_genome_wide[which(DTU_genome_wide$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
              
              if(Condition_DEBUG == 1)
              {
                cat("DTU_genome_wide_SIG_0\n")
                cat(str(DTU_genome_wide_SIG))
                cat("\n")
              }
              
              n_tested_GW<-length(unique(DTU_genome_wide$ensembl_gene_id))
              n_SIG_GW<-length(unique(DTU_genome_wide_SIG$ensembl_gene_id))
              Perc_SIG_GW<-round(100*(n_SIG_GW/n_tested_GW),2)
              
             
              
              if(dim(DTU_genome_wide_SIG)[1] >0)
              {
                DTU_genome_wide_SIG$n_tested_GW<-n_tested_GW
                DTU_genome_wide_SIG$n_SIG_GW<-n_SIG_GW
                DTU_genome_wide_SIG$Perc_SIG_GW<-Perc_SIG_GW
                
                if(Condition_DEBUG == 1)
                {
                  cat("DTU_genome_wide_SIG_1\n")
                  cat(str(DTU_genome_wide_SIG))
                  cat("\n")
                }
                
                DTU_genome_wide_SIG$VAR<-SELECTED_VARS_UPDATED_sel
                DTU_genome_wide_SIG$Proxy_VAR<-Proxy_array_sel
                
                
                if(Condition_DEBUG == 1)
                {
                  cat("DTU_genome_wide_SIG_2\n")
                  cat(str(DTU_genome_wide_SIG))
                  cat("\n")
                }
                #######################################
               
                list_GW[[k]]<-DTU_genome_wide_SIG
                
                # list_Block<-list()
                # list_CIS<-list()
                
              }#dim(DTU_genome_wide_SIG)[1] >0
            }#file.exists(DTU_genome_wide.rds)
            
            ### DTU block
            
            
            filename="DTU_Block.rds"
            
            n_SIG_Block<-0
            n_tested_Block<-0
            Perc_SIG_Block<-NA
            
            if (file.exists(filename)) {
              DTU_BLOCK = readRDS(filename)
              
              if(Condition_DEBUG == 1)
              {
                cat("DTU_BLOCK_INTIAL\n")
                cat(str(DTU_BLOCK))
                cat("\n")
              }
              
              DTU_BLOCK_SIG<-DTU_BLOCK[which(DTU_BLOCK$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
              
              if(Condition_DEBUG == 1)
              {
                cat("DTU_BLOCK_SIG_0\n")
                cat(str(DTU_BLOCK_SIG))
                cat("\n")
              }
              
              n_tested_Block<-length(unique(DTU_BLOCK$ensembl_gene_id))
              n_SIG_Block<-length(unique(DTU_BLOCK_SIG$ensembl_gene_id))
              Perc_SIG_Block<-round(100*(n_SIG_Block/n_tested_Block),2)
              
              DTU_BLOCK$n_tested_Block<-n_tested_Block
              DTU_BLOCK$n_SIG_Block<-n_SIG_Block
              DTU_BLOCK$Perc_SIG_Block<-Perc_SIG_Block
              
              
              # cat(str(DTU_BLOCK))
              # cat("\n")
              
              if(dim(DTU_BLOCK)[1] >0)
              {
               
                DTU_BLOCK$VAR<-SELECTED_VARS_UPDATED_sel
                DTU_BLOCK$Proxy_VAR<-Proxy_array_sel
                
                if(Condition_DEBUG == 1)
                {
                  cat("DTU_BLOCK_2\n")
                  cat(str(DTU_BLOCK))
                  cat("\n")
                }
                
                list_Block[[k]]<-DTU_BLOCK
                # quit(status=1)
                
                # list_CIS<-list()
                
              }
              
              
            }
            
            # cat("DTU_BLOCK_FINAL\n")
            
            
            ### CIS gene
            
            
            filename="DTU_CIS_gene.rds"
            
            n_SIG_CIS<-0
            n_tested_CIS<-0
            Perc_SIG_CIS<-NA
            
            if (file.exists(filename)) {
              DTU_CIS_gene = readRDS(filename)
              
              if(Condition_DEBUG == 1)
              {
                cat("DTU_CIS_gene\n")
                cat(str(DTU_CIS_gene))
                cat("\n")
              }
              
              DTU_CIS_gene_SIG<-DTU_CIS_gene[which(DTU_CIS_gene$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
              
              if(Condition_DEBUG == 1)
              {
                cat("DTU_CIS_gene_SIG_0\n")
                cat(str(DTU_CIS_gene_SIG))
                cat("\n")
              }
              
              n_tested_CIS<-length(unique(DTU_CIS_gene$ensembl_gene_id))
              n_SIG_CIS<-length(unique(DTU_CIS_gene_SIG$ensembl_gene_id))
              Perc_SIG_CIS<-round(100*(n_SIG_CIS/n_tested_CIS),2)
              
              DTU_CIS_gene$n_tested_CIS<-n_tested_CIS
              DTU_CIS_gene$n_SIG_CIS<-n_SIG_CIS
              DTU_CIS_gene$Perc_SIG_CIS<-Perc_SIG_CIS
              
              
              # cat(str(DTU_CIS_gene))
              # cat("\n")
              
              if(dim(DTU_CIS_gene)[1] >0)
              {
                DTU_CIS_gene$VAR<-SELECTED_VARS_UPDATED_sel
                DTU_CIS_gene$Proxy_VAR<-Proxy_array_sel
                
                if(Condition_DEBUG == 1)
                {
                  cat("DTU_CIS_gene_2\n")
                  cat(str(DTU_CIS_gene))
                  cat("\n")
                }
                
                list_CIS[[k]]<-DTU_CIS_gene
                
                
                # list_CIS<-list()
                
              }
              
              #  quit(status = 1)
              
              
            }
            
          }#file.exists(path10)
        }#k in 1:length(Proxy_array)
        

        if(length(list_CIS) >0)
        {
          
          Results_CIS = as.data.frame(data.table::rbindlist(list_CIS, fill=T), stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_CIS_0\n")
            cat(str(Results_CIS))
            cat("\n")
            cat(str(unique(Results_CIS$VAR)))
            cat("\n")
            cat(str(unique(Results_CIS$Proxy_VAR)))
            cat("\n")
          }
          
          list_CIS_DEF[[i]]<-Results_CIS
          
          
        }# length(list_CIS) >0
        
        if(length(list_Block) >0)
        {
          
          Results_Block = as.data.frame(data.table::rbindlist(list_Block, fill=T), stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_Block_0\n")
            cat(str(Results_Block))
            cat("\n")
            cat(str(unique(Results_Block$VAR)))
            cat("\n")
            cat(str(unique(Results_Block$Proxy_VAR)))
            cat("\n")
          }
          
          list_Block_DEF[[i]]<-Results_Block
          
        }# length(list_Block) >0
        
        if(length(list_GW) >0)
        {
          
          Results_GW = as.data.frame(data.table::rbindlist(list_GW, fill=T), stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_GW_0\n")
            cat(str(Results_GW))
            cat("\n")
            cat(str(unique(Results_GW$VAR)))
            cat("\n")
            cat(str(unique(Results_GW$Proxy_VAR)))
            cat("\n")
          }
          
          
          list_GW_DEF[[i]]<-Results_GW
        }# length(list_GW) >0
      }#dim(Proxy_file_UPDATED_sel)[1] >0
    }#i in 1:length(SELECTED_VARS_UPDATED)
    
    Condition_DEBUG <- 1
    
    Results_DEF<-data.frame()
    
    if(length(list_CIS_DEF) >0)
    {
      
      Results_CIS = as.data.frame(data.table::rbindlist(list_CIS_DEF, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_CIS_0\n")
        cat(str(Results_CIS))
        cat("\n")
        cat(str(unique(Results_CIS$VAR)))
        cat("\n")
      }
      

      Results_CIS_SIG<-Results_CIS[which(Results_CIS$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_CIS_SIG_0\n")
        cat(str(Results_CIS_SIG))
        cat("\n")
        cat(str(unique(Results_CIS_SIG$VAR)))
        cat("\n")
      }
      
      
      
      indx.int<-c(which(colnames(Results_CIS) == "VAR"),which(colnames(Results_CIS) == "Proxy_VAR"),which(colnames(Results_CIS) == "comparison"),
                  which(colnames(Results_CIS) == "ensembl_gene_id"),which(colnames(Results_CIS) == "HGNC"),which(colnames(Results_CIS) == "transcript_id"),
                  which(colnames(Results_CIS) == "ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS"),which(colnames(Results_CIS) == "coefficient_Haplotypes_specific_CELL_COUNTS"),
                  which(colnames(Results_CIS) == "n_breakdown_string"),
                  which(colnames(Results_CIS) == "n_SIG_CIS"),which(colnames(Results_CIS) == "n_tested_CIS"),which(colnames(Results_CIS) == "Perc_SIG_CIS"))
      
      Results_CIS_subset<-unique(Results_CIS[,indx.int])
      
      colnames(Results_CIS_subset)[which(colnames(Results_CIS_subset) == "ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS")]<-"CIS_gene_minuslogpvalue"
      colnames(Results_CIS_subset)[which(colnames(Results_CIS_subset) == "coefficient_Haplotypes_specific_CELL_COUNTS")]<-"CIS_Beta"
      colnames(Results_CIS_subset)[which(colnames(Results_CIS_subset) == "n_breakdown_string")]<-"CIS_n_breakdown_string"
      
      Results_DEF<-rbind(Results_CIS_subset,Results_DEF)
      
      
      
      
    }# length(list_CIS_DEF) >0
    
    if(length(list_Block_DEF) >0)
    {
      
      Results_Block = as.data.frame(data.table::rbindlist(list_Block_DEF, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_Block_0\n")
        cat(str(Results_Block))
        cat("\n")
        cat(str(unique(Results_Block$VAR)))
        cat("\n")
      }
      
      
      Results_Block_SIG<-Results_Block[which(Results_Block$ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS >= 1.3),]
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_Block_SIG_0\n")
        cat(str(Results_Block_SIG))
        cat("\n")
        cat(str(unique(Results_Block_SIG$VAR)))
        cat("\n")
      }
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_0\n")
        cat(str(Results_DEF))
        cat("\n")
        cat(str(unique(Results_DEF$VAR)))
        cat("\n")
      }
      
      
      
      indx.int<-c(which(colnames(Results_Block) == "VAR"),which(colnames(Results_Block) == "Proxy_VAR"),which(colnames(Results_Block) == "comparison"),
                  which(colnames(Results_Block) == "ensembl_gene_id"),which(colnames(Results_Block) == "HGNC"),which(colnames(Results_Block) == "transcript_id"),
                  which(colnames(Results_Block) == "ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS"),which(colnames(Results_Block) == "coefficient_Haplotypes_specific_CELL_COUNTS"),
                  which(colnames(Results_Block) == "n_breakdown_string"),
                  which(colnames(Results_Block) == "n_SIG_Block"),which(colnames(Results_Block) == "n_tested_Block"),which(colnames(Results_Block) == "Perc_SIG_Block"))
      
      Results_Block_subset<-unique(Results_Block[,indx.int])
      
      colnames(Results_Block_subset)[which(colnames(Results_Block_subset) == "ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS")]<-"Block_PCHiC_minuslogpvalue"
      colnames(Results_Block_subset)[which(colnames(Results_Block_subset) == "coefficient_Haplotypes_specific_CELL_COUNTS")]<-"Block_PCHiC_Beta"
      colnames(Results_Block_subset)[which(colnames(Results_Block_subset) == "n_breakdown_string")]<-"Block_PCHiC_n_breakdown_string"
      
      if(length(list_CIS_DEF) >0)
      {
        Results_DEF<-merge(Results_DEF,
                           Results_Block_subset,
                           by=c("VAR","Proxy_VAR","comparison","ensembl_gene_id","HGNC","transcript_id"),
                           all=T)
      }else{
        
        Results_DEF<-rbind(Results_Block_subset,Results_DEF)
        
      }
      
    }# length(list_Block_DEF) >0
    
    if(length(list_GW_DEF) >0)
    {
      
      Results_GW = as.data.frame(data.table::rbindlist(list_GW_DEF, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_GW_0\n")
        cat(str(Results_GW))
        cat("\n")
        cat(str(unique(Results_GW$VAR)))
        cat("\n")
      }
      
      
      
      indx.int<-c(which(colnames(Results_GW) == "VAR"),which(colnames(Results_GW) == "Proxy_VAR"),which(colnames(Results_GW) == "comparison"),
                  which(colnames(Results_GW) == "ensembl_gene_id"),which(colnames(Results_GW) == "HGNC"),which(colnames(Results_GW) == "transcript_id"),
                  which(colnames(Results_GW) == "ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS"),which(colnames(Results_GW) == "coefficient_Haplotypes_specific_CELL_COUNTS"),
                  which(colnames(Results_GW) == "n_breakdown_string"),
                  which(colnames(Results_GW) == "n_SIG_GW"),which(colnames(Results_GW) == "n_tested_GW"),which(colnames(Results_GW) == "Perc_SIG_GW"))
      
      Results_GW_subset<-unique(Results_GW[,indx.int])
      
      colnames(Results_GW_subset)[which(colnames(Results_GW_subset) == "ajusted.minuslogpvalue_Haplotypes_specific_CELL_COUNTS")]<-"Genome_wide_minuslogpvalue"
      colnames(Results_GW_subset)[which(colnames(Results_GW_subset) == "coefficient_Haplotypes_specific_CELL_COUNTS")]<-"Genome_wide_Beta"
      colnames(Results_GW_subset)[which(colnames(Results_GW_subset) == "n_breakdown_string")]<-"Genome_wide_n_breakdown_string"
      
      
      Results_GW_subset_SIG<-Results_GW_subset[which(Results_GW_subset$Genome_wide_minuslogpvalue >= 1.3),]
      
      Results_DEF<-merge(Results_DEF,
                         Results_GW_subset_SIG,
                         by=c("VAR","Proxy_VAR","comparison","ensembl_gene_id","HGNC","transcript_id"),
                         all=T)
      
      
      
      
      
    }# length(list_GW_DEF) >0
    
    if(dim(Results_DEF)[1] >0)
    {
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_BORN\n")
        cat(str(Results_DEF))
        cat("\n")
        cat(str(unique(Results_DEF$VAR)))
        cat("\n")
      }
      
      # ##################################################################################
      # quit(status = 1)
      
     
      
      path6<-paste(out,'FINAL_RESULTS','/', sep='')
      
      # cat("path6\n")
      # cat(sprintf(as.character(path6)))
      # cat("\n")
      
      
      if (file.exists(path6)){
        
        
        
        
      } else {
        dir.create(file.path(path6))
        
      }
      
      setwd(path6)
      
      ### Z-score the effect sizes ----
      
      Z_score_effect_sizes_BLOCK<-c(Results_DEF$CIS_Beta[!is.na(Results_DEF$CIS_Beta)],
                                    Results_DEF$Block_PCHiC_Beta[!is.na(Results_DEF$Block_PCHiC_Beta)])
      
      cat("Z_score_effect_sizes_BLOCK\n")
      cat(str(Z_score_effect_sizes_BLOCK))
      cat("\n")
      
      mean_Z_score_effect_sizes_BLOCK<-mean(Z_score_effect_sizes_BLOCK)
      
      cat("mean_Z_score_effect_sizes_BLOCK\n")
      cat(str(mean_Z_score_effect_sizes_BLOCK))
      cat("\n")
      
      sd_Z_score_effect_sizes_BLOCK<-sd(Z_score_effect_sizes_BLOCK)
      
      cat("sd_Z_score_effect_sizes_BLOCK\n")
      cat(str(sd_Z_score_effect_sizes_BLOCK))
      cat("\n")
      
      Results_DEF$CIS_Beta_Z_score<-(Results_DEF$CIS_Beta-mean_Z_score_effect_sizes_BLOCK)/sd_Z_score_effect_sizes_BLOCK
      Results_DEF$Block_PCHiC_Beta_Z_score<-(Results_DEF$Block_PCHiC_Beta-mean_Z_score_effect_sizes_BLOCK)/sd_Z_score_effect_sizes_BLOCK
      
      Z_score_effect_sizes_GW<-c(Results_DEF$Genome_wide_Beta[!is.na(Results_DEF$Genome_wide_Beta)])
      
      cat("Z_score_effect_sizes_GW\n")
      cat(str(Z_score_effect_sizes_GW))
      cat("\n")
      
      mean_Z_score_effect_sizes_GW<-mean(Z_score_effect_sizes_GW)
      
      cat("mean_Z_score_effect_sizes_GW\n")
      cat(str(mean_Z_score_effect_sizes_GW))
      cat("\n")
      
      sd_Z_score_effect_sizes_GW<-sd(Z_score_effect_sizes_GW)
      
      cat("sd_Z_score_effect_sizes_GW\n")
      cat(str(sd_Z_score_effect_sizes_GW))
      cat("\n")
      
      Results_DEF$Genome_wide_Beta_Z_score<-(Results_DEF$Genome_wide_Beta-mean_Z_score_effect_sizes_GW)/sd_Z_score_effect_sizes_GW
      
      
      cat("Results_DEF_FINAL_POST_Z\n")
      cat(str(Results_DEF))
      cat("\n")
      cat(str(unique(Results_DEF$VAR)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Results_DEF$CIS_Beta_Z_score)))))
      cat("\n")
      cat(sprintf(as.character(summary(Results_DEF$CIS_Beta_Z_score))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Results_DEF$Block_PCHiC_Beta_Z_score)))))
      cat("\n")
      cat(sprintf(as.character(summary(Results_DEF$Block_PCHiC_Beta_Z_score))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Results_DEF$Genome_wide_Beta_Z_score)))))
      cat("\n")
      cat(sprintf(as.character(summary(Results_DEF$Genome_wide_Beta_Z_score))))
      cat("\n")
      
      
      write.table(file="DTU_LM_logRatio_results_Haplotypes.tsv", Results_DEF, sep="\t", row.names = F,quote = F)

      
    }# length(list_DEF) >0
    
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
    make_option(c("--SELECTED_VARS_UPDATED "), type="character", default=NULL,
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
