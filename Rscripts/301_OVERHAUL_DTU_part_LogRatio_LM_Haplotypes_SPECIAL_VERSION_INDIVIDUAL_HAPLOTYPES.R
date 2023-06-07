

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


Filter_irrelevant_transcripts = function (option_list)
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
  
  #### Read transposed expression ----
  
  
  
  Transposed_Isoform_Expression<-readRDS(file=opt$Transposed_Isoform_Expression)
  colnames(Transposed_Isoform_Expression)[which(colnames(Transposed_Isoform_Expression) == "Sample_id")]<-"sample_id"
  
  # as.data.frame(fread(file=opt$Transposed_Isoform_Expression, sep=",", header=T) , stringsAsFactors=F)
  
  
  cat("Transposed_Isoform_Expression\n")
  cat(str(Transposed_Isoform_Expression))
  cat("\n")
  
  ### Read TRANSCRIPTS_table----
  
  
  TRANSCRIPTS_table<-as.data.frame(fread(file=opt$TRANSCRIPTS_table, sep="\t", header=F) , stringsAsFactors=F)
  
  colnames(TRANSCRIPTS_table)<-c("ensembl_gene_id","HGNC","transcript_id","transcript_name")
  
  cat("TRANSCRIPTS_table\n")
  cat(str(TRANSCRIPTS_table))
  cat("\n")
  
  ENST_array<-unique(colnames(Transposed_Isoform_Expression)[-which(colnames(Transposed_Isoform_Expression) == "sample_id")])
  
  cat("ENST_array\n")
  cat(str(ENST_array))
  cat("\n")
  
  ENSG_array<-unique(TRANSCRIPTS_table$ensembl_gene_id[which(TRANSCRIPTS_table$transcript_id%in%ENST_array)])
  
  cat("ENSG_array_0\n")
  cat(str(ENSG_array))
  cat("\n")
  
  # ENSG_array<-ENSG_array[which(ENSG_array == "ENSG00000000457")]
  # 
  # cat("ENSG_array_1\n")
  # cat(str(ENSG_array))
  # cat("\n")
  
  TRANSCRIPTS_table_FILTERED<-TRANSCRIPTS_table[which(TRANSCRIPTS_table$transcript_id%in%colnames(Transposed_Isoform_Expression)),]
  
  
  
  cat("TRANSCRIPTS_table_FILTERED_0\n")
  cat(str(TRANSCRIPTS_table_FILTERED))
  cat("\n")
  cat(str(unique(TRANSCRIPTS_table_FILTERED$transcript_id)))
  cat("\n")
  cat(str(unique(TRANSCRIPTS_table_FILTERED$ensembl_gene_id)))
  cat("\n")
  
  #### MASTER LOOP ----
  
  
  list_ACCEPTED_ENST<-list()
  
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
                
                
                
                Condition_DEBUG <- 0
                
                list_ACCEPTED_ENST<-list()
                
                for(i in 1:length(ENSG_array))
                {
                  ENSG_array_sel<-ENSG_array[i]
                  
                  cat("---------------->\t")
                  cat(sprintf(as.character(i)))
                  cat("\t")
                  cat(sprintf(as.character(ENSG_array_sel)))
                  cat("\n")
                  
                  TRANSCRIPTS_table_sel<-TRANSCRIPTS_table[which(TRANSCRIPTS_table$ensembl_gene_id == ENSG_array_sel),]
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("TRANSCRIPTS_table_sel\n")
                    cat(str(TRANSCRIPTS_table_sel))
                    cat("\n")
                  }
                  
                  Transposed_Isoform_Expression_sel<-Transposed_Isoform_Expression[,c(which(colnames(Transposed_Isoform_Expression) == "sample_id"),
                                                                                      which(colnames(Transposed_Isoform_Expression) %in% TRANSCRIPTS_table_sel$transcript_id))]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Transposed_Isoform_Expression_sel_0\n")
                    cat(str(Transposed_Isoform_Expression_sel))
                    cat("\n")
                  }
                  
                  
                  n_transcripts_per_gene<-dim(Transposed_Isoform_Expression_sel)[2]-1
                  
                  if(n_transcripts_per_gene > 1)
                  {
                    ################# melt the Transposed_Isoform_Expression_sel
                    
                    Transposed_Isoform_Expression_sel.m<-melt(Transposed_Isoform_Expression_sel, id.vars=c("sample_id"),value.name = "log2TPM",variable.name = "transcript_id")
                    Transposed_Isoform_Expression_sel.m$TPM<- -1+2^(Transposed_Isoform_Expression_sel.m$log2TPM)
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_0\n")
                      cat(str(Transposed_Isoform_Expression_sel.m))
                      cat("\n")
                      
                      cat("sample_id\n")
                      cat(str(unique(Transposed_Isoform_Expression_sel.m$sample_id)))
                      cat("\n")
                      
                    }
                    
                    
                    Transposed_Isoform_Expression_sel.m_ADAPTED<-merge(Transposed_Isoform_Expression_sel.m,
                                                                       Haplotype_PEER_G_HET_haplotypes_sel,
                                                                       by="sample_id")
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_ADAPTED_0\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_ADAPTED))
                      cat("\n")
                      
                      cat("sample_id\n")
                      cat(str(unique(Transposed_Isoform_Expression_sel.m_ADAPTED$sample_id)))
                      cat("\n")
                      
                      cat(sprintf(as.character(names(summary(Transposed_Isoform_Expression_sel.m_ADAPTED$Haplotype)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Transposed_Isoform_Expression_sel.m_ADAPTED$Haplotype))))
                      cat("\n")
                      
                    }
                    
                    
                    Transposed_Isoform_Expression_sel.m_ADAPTED.dt<-data.table(Transposed_Isoform_Expression_sel.m_ADAPTED,
                                                                               key=c("Haplotype","transcript_id"))
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_ADAPTED.dt_2\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_ADAPTED.dt))
                      cat("\n")
                    }
                    
                    ### calculate median TPM per Haplotype
                    
                    Summary_table<-as.data.frame(Transposed_Isoform_Expression_sel.m_ADAPTED.dt[, .(median=round(as.numeric(summary(TPM)[3]),3)),
                                                                                                by=key(Transposed_Isoform_Expression_sel.m_ADAPTED.dt)],stringsAsFactors=F)
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Summary_table_0\n")
                      cat(str(Summary_table))
                      cat("\n")
                      
                      #  quit(status = 1)
                    }
                    
                    Summary_table.dt<-data.table(Summary_table,
                                                 key=c("Haplotype"))
                    
                    #### calculate total GENE EXP per Haplotype as the sum of the medians
                    
                    Summary_table_TOTAL<-as.data.frame(Summary_table.dt[, .(TOTAL_GENE_EXP_median=sum(median)),
                                                                        by=key(Summary_table.dt)],stringsAsFactors=F)
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Summary_table_TOTAL_0\n")
                      cat(str(Summary_table_TOTAL))
                      cat("\n")
                      
                      #  quit(status = 1)
                    }
                    
                    Summary_table<-merge(Summary_table,
                                         Summary_table_TOTAL,
                                         by=c("Haplotype"))
                    
                    #### Calculate Transcript Ratio
                    
                    Summary_table$Transcript_Ratio<-(Summary_table$median/Summary_table$TOTAL_GENE_EXP_median)
                    
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Summary_table_1\n")
                      cat(str(Summary_table))
                      cat("\n")
                      
                      # quit(status = 1)
                    }
                    
                    
                    Summary_table.dt<-data.table(Summary_table,
                                                 key=c("transcript_id"))
                    
                    #### Obtain the max transcript ratio per transcript from all the Haplotypes
                    
                    
                    Summary_table_FILTERED<-as.data.frame(Summary_table.dt[,.SD[which.max(Transcript_Ratio)],
                                                                           by=key(Summary_table.dt)],stringsAsFactors=F)
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Summary_table_FILTERED_0\n")
                      cat(str(Summary_table_FILTERED))
                      cat("\n")
                      
                    }
                    
                    ### Filter if it is less than 10%
                    
                    if(SELECTED_VARS_UPDATED_sel == "chr9_135920196_C_T")
                    {
                      Summary_table_FILTERED<-Summary_table_FILTERED[which(Summary_table_FILTERED$Transcript_Ratio >= 0.01),]
                      
                    }else{
                      
                      Summary_table_FILTERED<-Summary_table_FILTERED[which(Summary_table_FILTERED$Transcript_Ratio >= 0.1),]
                      
                    }
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Summary_table_FILTERED_1\n")
                      cat(str(Summary_table_FILTERED))
                      cat("\n")
                    }
                    
                    if(dim(Summary_table_FILTERED)[1] >1)
                    {
                      # INTERVAL_isoform_EXP_Filtered_sel_restricted.m_FILTERED<-INTERVAL_isoform_EXP_Filtered_sel_restricted.m[which(INTERVAL_isoform_EXP_Filtered_sel_restricted.m$transcript_id%in%Summary_table_FILTERED$transcript_id),]
                      
                      PASSED_ENST<-unique(as.character(Summary_table_FILTERED$transcript_id))
                      
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("PASSED_ENST_0\n")
                        cat(str(PASSED_ENST))
                        cat("\n")
                        
                        # ###################################################################################################################
                        # quit(status = 1)
                      }
                      
                      list_ACCEPTED_ENST[[i]]<-PASSED_ENST
                      
                      
                    }else{
                      
                      # After filtering the TR > 10% int at least 1 Haplotype the gene only has 1 isoform
                    }
                  }#n_transcripts_per_gene > 1
                }#i ENSG_array   
                
                Condition_DEBUG <- 1
                
                
                if(length(list_ACCEPTED_ENST) >0)
                {
                  
                  FINAL_transcripts = unique(unlist(list_ACCEPTED_ENST))
                  
                  FINAL_transcripts<-FINAL_transcripts[!is.null(FINAL_transcripts)]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("FINAL_transcripts_1\n")
                    cat(str(FINAL_transcripts))
                    cat("\n")
                  }
                  
                  TRANSCRIPTS_table_FILTERED_FINAL<-unique(TRANSCRIPTS_table_FILTERED[which(TRANSCRIPTS_table_FILTERED$transcript_id%in%FINAL_transcripts),])
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("TRANSCRIPTS_table_FILTERED_FINAL:1\n")
                    cat(str(TRANSCRIPTS_table_FILTERED_FINAL))
                    cat("\n")
                  }
                  
                  
                  path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                  
                  setwd(path10)
                  
                  
                  saveRDS(TRANSCRIPTS_table_FILTERED_FINAL,file=paste("DTU_PASS_Transcripts_",comparison_2,'.rds', sep=''))
                  
                  
                }# length(list_ACCEPTED_ENST) >0
              }# l in 1:length(array_haplotypes)
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



LogRatio_LM_model = function (option_list)
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
  
  #### Read transposed expression ----
  
  
  
  Transposed_Isoform_Expression<-readRDS(file=opt$Transposed_Isoform_Expression)
  colnames(Transposed_Isoform_Expression)[which(colnames(Transposed_Isoform_Expression) == "Sample_id")]<-"sample_id"
  
  # as.data.frame(fread(file=opt$Transposed_Isoform_Expression, sep=",", header=T) , stringsAsFactors=F)
  
  
  cat("Transposed_Isoform_Expression\n")
  cat(str(Transposed_Isoform_Expression))
  cat("\n")
  
  ### Read TRANSCRIPTS_table----
  
  
  TRANSCRIPTS_table<-as.data.frame(fread(file=opt$TRANSCRIPTS_table, sep="\t", header=F) , stringsAsFactors=F)
  
  colnames(TRANSCRIPTS_table)<-c("ensembl_gene_id","HGNC","transcript_id","transcript_name")
  
  cat("TRANSCRIPTS_table\n")
  cat(str(TRANSCRIPTS_table))
  cat("\n")
  
  ENST_array<-unique(colnames(Transposed_Isoform_Expression)[-which(colnames(Transposed_Isoform_Expression) == "sample_id")])
  
  cat("ENST_array\n")
  cat(str(ENST_array))
  cat("\n")
  
  ENSG_array<-unique(TRANSCRIPTS_table$ensembl_gene_id[which(TRANSCRIPTS_table$transcript_id%in%ENST_array)])
  
  cat("ENSG_array_0\n")
  cat(str(ENSG_array))
  cat("\n")
  
  
  
  TRANSCRIPTS_table_FILTERED<-TRANSCRIPTS_table[which(TRANSCRIPTS_table$transcript_id%in%colnames(Transposed_Isoform_Expression)),]
  
  
  
  cat("TRANSCRIPTS_table_FILTERED_0\n")
  cat(str(TRANSCRIPTS_table_FILTERED))
  cat("\n")
  cat(str(unique(TRANSCRIPTS_table_FILTERED$transcript_id)))
  cat("\n")
  cat(str(unique(TRANSCRIPTS_table_FILTERED$ensembl_gene_id)))
  cat("\n")
  
  ############################ MASTER LOOP ##########################################################
  
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
      
      Condition_DEBUG <- 1
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
                
                #### DTU_PASS_Transcripts
                
                setwd(path10)
                
                TRANSCRIPTS_table_FILTERED_FINAL<-readRDS(file=paste("DTU_PASS_Transcripts_",comparison_2,'.rds', sep=''))
                
                if(Condition_DEBUG == 1)
                {
                  cat("TRANSCRIPTS_table_FILTERED_FINAL\n")
                  cat(str(TRANSCRIPTS_table_FILTERED_FINAL))
                  cat("\n")
                }
                
                ### Residuals file
                
                setwd(path10)
                
                filename<-paste("DTU_LM_HET_RESIDUALS_",comparison_2,'.csv',sep='')
                
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
                
                
                
                ENSG_array<-unique(TRANSCRIPTS_table_FILTERED_FINAL$ensembl_gene_id)
                
               
                
                # ENSG_array <- "ENSG00000000457"
                
                
                # if(Condition_DEBUG == 1)
                # {
                cat("ENSG_array\n")
                cat(str(ENSG_array))
                cat("\n")
                # }
                
                
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
                  
                  
                  
                  
                  TRANSCRIPTS_table_FILTERED_FINAL_sel<-TRANSCRIPTS_table_FILTERED_FINAL[which(TRANSCRIPTS_table_FILTERED_FINAL$ensembl_gene_id%in%ENSG_array_sel),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("TRANSCRIPTS_table_FILTERED_FINAL_sel_0\n")
                    cat(str(TRANSCRIPTS_table_FILTERED_FINAL_sel))
                    cat("\n")
                  }
                  
                  
                  Transposed_Isoform_Expression_sel<-Transposed_Isoform_Expression[,c(which(colnames(Transposed_Isoform_Expression) == "sample_id"),
                                                                                      which(colnames(Transposed_Isoform_Expression) %in% TRANSCRIPTS_table_FILTERED_FINAL_sel$transcript_id))]
                  if(Condition_DEBUG == 1)
                  {
                    cat("Transposed_Isoform_Expression_sel_1\n")
                    cat(str(Transposed_Isoform_Expression_sel))
                    cat("\n")
                  }
                  
                  ##### MASTER merge -----
                  
                  
                  
                  
                  
                  n_transcripts_per_gene<-dim(Transposed_Isoform_Expression_sel)[2]-1
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("n_transcripts_per_gene\n")
                    cat(str(n_transcripts_per_gene))
                    cat("\n")
                  }
                  
                  if(n_transcripts_per_gene > 1)
                  {
                    ################# melt the Transposed_Isoform_Expression_sel
                    
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
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED<-Transposed_Isoform_Expression_sel.m[which(Transposed_Isoform_Expression_sel.m$sample_id%in%Haplotype_PEER_G_HET_haplotypes_sel$sample_id),]
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED_0\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED))
                      cat("\n")
                    }
                    
                    
                    #### calculate log model ----
                    
                    #### Find lowest expression per transcript that is different from 0
                    
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO<-Transposed_Isoform_Expression_sel.m_REDUCED[which(Transposed_Isoform_Expression_sel.m_REDUCED$TPM >0),]
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO_0\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO))
                      cat("\n")
                    }
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED_ZERO<-droplevels(Transposed_Isoform_Expression_sel.m_REDUCED[which(Transposed_Isoform_Expression_sel.m_REDUCED$TPM <= 0),])
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED_ZERO_0\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO))
                      cat("\n")
                      cat(sprintf(as.character(names(summary(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$TPM)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$TPM))))
                      cat("\n")
                    }
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO,
                                                                                       key=c("transcript_id"))
                    
                    
                    
                    Zero_imputation_values<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO.dt[,.(min_TPM=min(TPM)),by=key(Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO.dt)], stringsAsFactors=F)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Zero_imputation_values_0\n")
                      cat(str(Zero_imputation_values))
                      cat("\n")
                      
                      
                    }
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED_ZERO<-merge(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO,
                                                                            Zero_imputation_values,
                                                                            by="transcript_id",
                                                                            all.x=T)
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$TPM<-0.65*Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$min_TPM
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED_ZERO_1\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO))
                      cat("\n")
                      cat(sprintf(as.character(names(summary(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$TPM)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO$TPM))))
                      cat("\n")
                      
                      # ##################################################
                      # quit(status = 1)
                    }
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED_ZERO<-Transposed_Isoform_Expression_sel.m_REDUCED_ZERO[,-which(colnames(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO)== "min_TPM")]
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED<-rbind(Transposed_Isoform_Expression_sel.m_REDUCED_ZERO,
                                                                       Transposed_Isoform_Expression_sel.m_REDUCED_NO_ZERO)
                    
                    
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED_RBIND\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED))
                      cat("\n")
                      cat(sprintf(as.character(names(summary(Transposed_Isoform_Expression_sel.m_REDUCED$TPM)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Transposed_Isoform_Expression_sel.m_REDUCED$TPM))))
                      cat("\n")
                    }
                    
                    
                    ### Calculate transcript_reference
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED,
                                                                               key=c("transcript_id"))
                    
                    Reference_value<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED.dt[,.(mean_TPM=mean(TPM)),by=key(Transposed_Isoform_Expression_sel.m_REDUCED.dt)], stringsAsFactors=F)
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Reference_value_0\n")
                      cat(str(Reference_value))
                      cat("\n")
                      
                      
                    }
                    
                    max_value<-max(Reference_value$mean_TPM)
                    
                    Reference_value_MAX<-Reference_value[which(Reference_value$mean_TPM == max_value),]
                    
                    Reference_transcript<-Reference_value_MAX$transcript_id[1]
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("max_value\n")
                      cat(sprintf(as.character(max_value)))
                      cat("\n")
                      
                      cat("Reference_value_MAX_0\n")
                      cat(str(Reference_value_MAX))
                      cat("\n")
                      
                      cat("Reference_transcript\n")
                      cat(sprintf(as.character(Reference_transcript)))
                      cat("\n")
                      
                      
                    }
                    
                    
                    ### calculate summatory TPM per sample
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED,
                                                                               key=c("sample_id"))
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED.dt_0\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED.dt))
                      cat("\n")
                    }
                    
                    
                    Summary_table_GENE_EXP<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED.dt[, .(sum_GENE_EXP=sum(TPM)),
                                                                                                         by=key(Transposed_Isoform_Expression_sel.m_REDUCED.dt)],stringsAsFactors=F)
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Summary_table_GENE_EXP_0\n")
                      cat(str(Summary_table_GENE_EXP))
                      cat("\n")
                      
                      #  quit(status = 1)
                    }
                    
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED<-merge(Transposed_Isoform_Expression_sel.m_REDUCED,
                                                                       Summary_table_GENE_EXP,
                                                                       by="sample_id")
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      cat("Transposed_Isoform_Expression_sel.m_REDUCED_1\n")
                      cat(str(Transposed_Isoform_Expression_sel.m_REDUCED))
                      cat("\n")
                    }
                    
                    Transposed_Isoform_Expression_sel.m_REDUCED.dt<-data.table(Transposed_Isoform_Expression_sel.m_REDUCED,
                                                                               key=c("sample_id","transcript_id"))
                    
                    Ratio_df<-as.data.frame(Transposed_Isoform_Expression_sel.m_REDUCED.dt[, .(TPM=TPM,
                                                                                               sum_GENE_EXP=sum_GENE_EXP,
                                                                                               Ratio=TPM/sum_GENE_EXP),by=key(Transposed_Isoform_Expression_sel.m_REDUCED.dt)],stringsAsFactors=F)
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Ratio_df_0\n")
                      cat(str(Ratio_df))
                      cat("\n")
                      
                      cat("distrib_ratios\n")
                      cat(sprintf(as.character(names(summary(Ratio_df$Ratio)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Ratio_df$Ratio))))
                      cat("\n")
                      
                    }
                    
                    
                    
                    
                    #### scale ratios to reference transcripts
                    
                    Ratio_df_reference<-droplevels(Ratio_df[which(Ratio_df$transcript_id == Reference_transcript),c(which(colnames(Ratio_df) == "sample_id"),
                                                                                                                    which(colnames(Ratio_df) == "transcript_id"),
                                                                                                                    which(colnames(Ratio_df) == "Ratio"))])
                    
                    colnames(Ratio_df_reference)[which(colnames(Ratio_df_reference) == "Ratio")]<-"Reference_ratio_value"
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Ratio_df_reference_0\n")
                      cat(str(Ratio_df_reference))
                      cat("\n")
                      
                      cat("distrib_ratios\n")
                      cat(sprintf(as.character(names(summary(Ratio_df_reference$Reference_ratio_value)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Ratio_df_reference$Reference_ratio_value))))
                      cat("\n")
                      
                    }
                    
                    Ratio_df_reference<-unique(Ratio_df_reference[,-which(colnames(Ratio_df_reference) == "transcript_id")])
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Ratio_df_reference_1\n")
                      cat(str(Ratio_df_reference))
                      cat("\n")
                    }
                    
                    Ratio_df<-merge(Ratio_df,
                                    Ratio_df_reference,
                                    by="sample_id",
                                    all.x=T)
                    
                    Ratio_df$scaled_ratio<-Ratio_df$Ratio/Ratio_df$Reference_ratio_value
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Ratio_df_1\n")
                      cat(str(Ratio_df))
                      cat("\n")
                      
                      cat("distrib_ratios\n")
                      cat(sprintf(as.character(names(summary(Ratio_df$scaled_ratio)))))
                      cat("\n")
                      cat(sprintf(as.character(summary(Ratio_df$scaled_ratio))))
                      cat("\n")
                      
                    }
                    
                    check<-Ratio_df[which(Ratio_df$scaled_ratio > 1),]
                    
                    if(dim(check)[1] >0)
                    {
                      check2<-Ratio_df[which(Ratio_df$sample_id%in%check$sample_id),]
                      check2<-check2[order(check2$sample_id),]
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("check2\n")
                        cat(str(check2))
                        cat("\n")
                      }
                      
                    }
                    
                    
                    
                    
                    #### calculate the log
                    
                    Ratio_df$logscaled_ratio<-log(Ratio_df$scaled_ratio)
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("Ratio_df_DEF\n")
                      cat(str(Ratio_df))
                      cat("\n")
                      
                      
                      
                    }
                    
                    ##### Fit the LM per transcript ----
                    
                    ENST_array<-unique(as.character(Ratio_df$transcript_id))
                    
                    
                    if(Condition_DEBUG == 1)
                    {
                      
                      cat("ENST_array_1\n")
                      cat(str(ENST_array))
                      cat("\n")
                      
                    }
                    
                    
                    list_transcript<-list()
                    
                    Condition_DEBUG <- 0
                    
                    
                    for(z in 1:length(ENST_array))
                    {
                      ENST_array_sel<-ENST_array[z]
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("--->\t")
                        cat(sprintf(as.character(ENST_array_sel)))
                        cat("\n")
                        
                      }
                      
                      Ratio_df_ENST_sel<-Ratio_df[which(Ratio_df$transcript_id%in%ENST_array_sel),] 
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("Ratio_df_ENST_sel_0\n")
                        cat(str(Ratio_df_ENST_sel))
                        cat("\n")
                        
                        
                        # quit(status = 1)
                        
                      }
                      
                      
                      
                      
                      
                      Ratio_df_ENST_sel<-merge(Ratio_df_ENST_sel,
                                               Haplotype_PEER_G_HET_haplotypes_sel,
                                               by="sample_id")
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("Ratio_df_ENST_sel_1\n")
                        cat(str(Ratio_df_ENST_sel))
                        cat("\n")
                        cat(sprintf(as.character(colnames(Ratio_df_ENST_sel))))
                        cat("\n")
                        
                        
                        # quit(status = 1)
                        
                      }
                      
                      n_summary<-as.numeric(summary(Ratio_df_ENST_sel$Haplotype))
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("n_summary_0\n")
                        cat(str(n_summary))
                        cat("\n")
                        
                        
                        # quit(status = 1)
                        
                      }
                      
                      n_summary_names<-names(summary(Ratio_df_ENST_sel$Haplotype))
                      
                      
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
                      ### FULL linear model of the logRatio (Haplotype+covs) ----
                      
                      unselected_columns<-c("sample_id","transcript_id","TPM","sum_GENE_EXP","Ratio","Reference_ratio_value","scaled_ratio")
                      
                      
                      
                      # Selected_columns<-c("logRatio",
                      #                     "Haplotype",
                      #                     "age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth",
                      #                     "sex",
                      #                     "Season_Winter","Season_Autumn","Season_Spring","Season_Summer",
                      #                     "sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5","sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                      #                     "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                      #                     "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24","PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35",
                      #                     "HCT_PCT___RNA","HGB_g_dL___RNA","IRF_PCT___RNA","MCH_pg___RNA","MCHC_g_dL___RNA","MCV_fL___RNA","RBC_10_12_L___RNA","RDW_SD_fL___RNA","RET_PCT___RNA")
                      
                      Ratio_df_ENST_sel_LM<-Ratio_df_ENST_sel[,-which(colnames(Ratio_df_ENST_sel)%in%unselected_columns)]
                      
                      row.names(Ratio_df_ENST_sel_LM)<-Ratio_df_ENST_sel$sample_id
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("Ratio_df_ENST_sel_LM_0\n")
                        cat(str(Ratio_df_ENST_sel_LM))
                        cat("\n")
                        
                        
                        
                        
                      }
                      
                      
                      model<-lm(logscaled_ratio ~ ., data=Ratio_df_ENST_sel_LM)
                      
                      
                      
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
                      
                      A.df<-as.data.frame(cbind(ENST_array_sel,pvalue_Haplotypes_1, coefficient_Haplotypes_1,RUN_OUT_OF_NAMES_DEF))
                      
                      colnames(A.df)<-c("transcript_id","pvalue_Haplotypes_specific_CELL_COUNTS","coefficient_Haplotypes_specific_CELL_COUNTS","n_breakdown_string")
                      
                      A.df$pvalue_Haplotypes_specific_CELL_COUNTS<-as.numeric(A.df$pvalue_Haplotypes_specific_CELL_COUNTS)
                      A.df$coefficient_Haplotypes_specific_CELL_COUNTS<-as.numeric(A.df$coefficient_Haplotypes_specific_CELL_COUNTS)
                      
                      list_transcript[[z]]<-A.df
                      
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("A.df\n")
                        cat(str(A.df))
                        cat("\n")
                        
                        
                        # quit(status=1)
                      }
                      
                      ### REDUCED linear model of the logRatio (covs) ----
                      
                      
                      
                      unselected_columns<-c("sample_id","transcript_id","sum_GENE_EXP","TPM","Reference_ratio_value","scaled_ratio","Haplotype","logscaled_ratio")
                      
                      
                      
                      # Selected_columns<-c("logRatio",
                      #                     "age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth",
                      #                     "sex",
                      #                     "Season_Winter","Season_Autumn","Season_Spring","Season_Summer",
                      #                     "sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5","sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                      #                     "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                      #                     "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24","PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35",
                      #                     "HCT_PCT___RNA","HGB_g_dL___RNA","IRF_PCT___RNA","MCH_pg___RNA","MCHC_g_dL___RNA","MCV_fL___RNA","RBC_10_12_L___RNA","RDW_SD_fL___RNA","RET_PCT___RNA")
                      
                      Ratio_df_ENST_sel_LM<-Ratio_df_ENST_sel[,-which(colnames(Ratio_df_ENST_sel)%in%unselected_columns)]
                      
                      row.names(Ratio_df_ENST_sel_LM)<-Ratio_df_ENST_sel$sample_id
                      
                      if(Condition_DEBUG == 1)
                      {
                        
                        cat("Ratio_df_ENST_sel_LM_REDUCED_LM\n")
                        cat(str(Ratio_df_ENST_sel_LM))
                        cat("\n")
                        
                      }
                      
                      
                      model<-lm(Ratio ~ ., data=Ratio_df_ENST_sel_LM)
                      
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
                        
                        
                       
                        
                        
                        residual_df$transcript_id<-ENST_array_sel
                        
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
                                                                    id_cols=c("transcript_id"),
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
                        
                        residual_df$transcript_id<-ENST_array_sel
                        
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
                                                                    id_cols=c("transcript_id"),
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
                      
                      
                      
                      
                      
                    }# z ENST_array
                    
                    Condition_DEBUG <- 0
                    
                    if(length(list_transcript) >0)
                    {
                      
                      Results_per_gene = as.data.frame(data.table::rbindlist(list_transcript, fill=T), stringsAsFactors=F)
                      
                      Results_per_gene$ensembl_gene_id<-ENSG_array_sel
                      
                      list_RESULT[[i]]<-Results_per_gene
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("Results_per_gene_0\n")
                        cat(str(Results_per_gene))
                        cat("\n")
                        
                        # ###################################################################
                        # quit(status = 1)
                        
                        #quit(status = 1)
                      }
                    }# length(list_transcript) >0
                    
                  }#n_transcripts_per_gene > 1
                  
                }#i ENSG_array
                
                Condition_DEBUG <- 1
                
                
                if(length(list_RESULT) >0)
                {
                  
                  Results_uncorrected = as.data.frame(data.table::rbindlist(list_RESULT, fill=T), stringsAsFactors=F)
                  
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Results_uncorrected_0\n")
                    cat(str(Results_uncorrected))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$ensembl_gene_id)))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$transcript_id)))
                    cat("\n")
                  }
                  
                  Results_uncorrected<-merge(Results_uncorrected,
                                             TRANSCRIPTS_table_FILTERED_FINAL,
                                             by=c("ensembl_gene_id","transcript_id"))
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Results_uncorrected_1\n")
                    cat(str(Results_uncorrected))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$ensembl_gene_id)))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$transcript_id)))
                    cat("\n")
                  }
                  
                  
                  Results_uncorrected$comparison<-comparison
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("Results_uncorrected_2\n")
                    cat(str(Results_uncorrected))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$ensembl_gene_id)))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$transcript_id)))
                    cat("\n")
                    cat(str(unique(Results_uncorrected$comparison)))
                    cat("\n")
                  }
                  
                  list_haplotypes[[l]]<-Results_uncorrected
                  
                }# length(list_RESULT) >0
              }# l in 1:length(array_haplotypes)
              
              Condition_DEBUG <- 1
              
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
                
                
                saveRDS(Results_uncorrected_all_comparisons,file=paste("DTU_LogLM_Haplotypes_RESULTS_NOMINAL_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'.rds', sep=''))
                
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
    make_option(c("--TRANSCRIPTS_table"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--INTERVAL_isoform_EXP"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Transposed_Isoform_Expression"), type="character", default=NULL,
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
  
  # Filter_irrelevant_transcripts(opt)
  LogRatio_LM_model(opt)
  
  
  
  
}


###########################################################################

system.time( main() )
