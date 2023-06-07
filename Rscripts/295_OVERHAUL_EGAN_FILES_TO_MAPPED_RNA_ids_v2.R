

suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
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
library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

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


opt = NULL


Data_wrangling = function(option_list)
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
  
  #### Read INTERVAL_covariates_and_PEER_factors ----
  
  INTERVAL_covariates_and_PEER_factors<-as.data.frame(fread(file=opt$INTERVAL_covariates_and_PEER_factors, sep="\t", header=T) , stringsAsFactors=F)
  
  
  
  # cat("INTERVAL_covariates_and_PEER_factors\n")
  # cat(str(INTERVAL_covariates_and_PEER_factors))
  # cat("\n")
  
  ### Read covariate_phenotype_correspondence ----
  
  
  covariate_phenotype_correspondence<-as.data.frame(fread(file=opt$covariate_phenotype_correspondence, sep="\t", header=T) , stringsAsFactors=F)
  
  
  cat("covariate_phenotype_correspondence\n")
  cat(str(covariate_phenotype_correspondence))
  cat("\n")
  
  
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  
  ALL_dB_Thresholded<-ALL_dB[which(ALL_dB$finemap_prob>= 0.1),]
  
  rm(ALL_dB)
  
  
  
  cat("ALL_dB_Thresholded_0\n")
  cat(str(ALL_dB_Thresholded))
  cat("\n")
  
 
  #### Selected_vars ----
  
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS))
  cat("\n")
  
  #### READ and transform Thomas_equivalence ----
  
  Thomas_equivalence<-as.data.frame(fread(opt$Thomas_equivalence,
                                   sep="\t", header=T), stringsAsFactors = F)
  
  
  cat("Thomas_equivalence_0\n")
  cat(str(Thomas_equivalence))
  cat("\n")
  cat(str(unique(Thomas_equivalence$sample_id)))
  cat("\n")
  
  #Thomas_equivalence<-Thomas_equivalence[-which(Thomas_equivalence$egan_id != "")]
  
 #### Read previous files ----
  
  setwd(out)
  
  Proxy_file<-as.data.frame(fread(file="Proxy_R2_TABLE_ALL.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("Proxy_file\n")
  cat(str(Proxy_file))
  cat("\n")
  cat(str(unique(Proxy_file$VAR)))
  cat("\n")
  
  
  VAR_and_Proxys<-as.data.frame(fread(file="ALL_variants.tsv",sep="\t",header=T), stringsAsFactors=F)
  
  cat("VAR_and_Proxys\n")
  cat(str(VAR_and_Proxys))
  cat("\n")
 
  # quit(status = 1)
  #### READ and transform EGAN_FILES
  
  EGAN_files_path<-paste(out,'EGAN_files','/',sep='')
  
  file_list <- list.files(path=EGAN_files_path, include.dirs = FALSE)
  indexes_sel <- grep("_EGAN_file\\.tsv$",file_list)
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  
  colnames(file_list_sel)<-"file"
  
  cat("file_list_sel\n")
  cat(str(file_list_sel))
  cat("\n")
  
  # quit(status = 1)
  
  EGAN_df<-data.frame()
  
  
  
  ABSENT_RNA_ids<-list()
  
  Condition_DEBUG <- 1
  
  for(i in 1:dim(file_list_sel)[1])
  {
    file_sel<-file_list_sel$file[i]
    
    cat("------------------->\t")
    cat(sprintf(as.character(file_sel)))
    cat("\n")
    
    setwd(EGAN_files_path)
    
    if(file.exists(file_sel))
    {
      
      tmp_vcf<-readLines(file_sel)
      
      cat("tmp_vcf:\n")
      cat(str(tmp_vcf))
      cat("\n")
      
      n_lines<-length(tmp_vcf)
      
      cat("n_lines:\n")
      cat(str(n_lines))
      cat("\n")
      
      indx_COMMENT<-grep("^##",tmp_vcf)
      
      # cat("indx_COMMENT:\n")
      # cat(str(indx_COMMENT))
      # cat("\n")
        
      indx_header<-grep("^#CHROM",tmp_vcf)
      
      cat("indx_header:\n")
      cat(str(indx_header))
      cat("\n")
      
      indx_lines<-grep("^[^#]",tmp_vcf)
      
      # cat("indx_lines:\n")
      # cat(str(indx_lines))
      # cat("\n")
      
      tmp_vcf_header<-tmp_vcf[indx_header]
      
      # cat("tmp_vcf_header_0:\n")
      # cat(str(tmp_vcf_header))
      # cat("\n")
      
      tmp_vcf_header<-gsub("^#","",tmp_vcf_header)
      
      # cat("tmp_vcf_header_1:\n")
      # cat(str(tmp_vcf_header))
      # cat("\n")
      
      tmp_vcf_header_terms<-unlist(strsplit(tmp_vcf_header, split="\t"))
      
      if(Condition_DEBUG == 1)
      {
        cat("tmp_vcf_header_terms_0:\n")
        cat(str(tmp_vcf_header_terms))
        cat("\n")
      }
      
      if(n_lines > indx_header)
      {
        EGAN_FILE<-as.data.frame(fread(file=file_sel, header=F, sep="\t", skip=indx_header), stringsAsFactors=F)
        
        colnames(EGAN_FILE)<-tmp_vcf_header_terms
        
        # cat("EGAN_FILE:\n")
        # cat(str(EGAN_FILE))
        # cat("\n")
        
        if(dim(EGAN_FILE)[1] >0)
        {
          EGAN_FILE$VAR_38<-paste(EGAN_FILE$CHROM,EGAN_FILE$POS,EGAN_FILE$REF,EGAN_FILE$ALT, sep="_")
          
          indx.dep<-c(which(colnames(EGAN_FILE) == "CHROM"),which(colnames(EGAN_FILE) == "POS"),which(colnames(EGAN_FILE) == "ID"),which(colnames(EGAN_FILE) == "REF"),which(colnames(EGAN_FILE) == "ALT"),
                      which(colnames(EGAN_FILE) == "QUAL"),which(colnames(EGAN_FILE) == "FILTER"),which(colnames(EGAN_FILE) == "INFO"),which(colnames(EGAN_FILE) == "FORMAT"))
          
          
          EGAN_FILE_subset<-unique(EGAN_FILE[,-indx.dep])
          
          # cat("EGAN_FILE_subset_0:\n")
          # cat(str(EGAN_FILE_subset))
          # cat("\n")
          # cat(str(unique(EGAN_FILE_subset$VAR_38)))
          # cat("\n")
          
          indx.VAR<-which(colnames(EGAN_FILE_subset) == "VAR_38")
          indx.no.VAR<-which(colnames(EGAN_FILE_subset) != "VAR_38")
          
          EGAN_FILE_subset<-EGAN_FILE_subset[,c(indx.VAR,indx.no.VAR)]
          
          # cat("EGAN_FILE_subset_1:\n")
          # cat(str(EGAN_FILE_subset))
          # cat("\n")
          # cat(str(unique(EGAN_FILE_subset$VAR_38)))
          # cat("\n")
          
          EGAN_FILE_subset.m<-reshape2::melt(EGAN_FILE_subset, id.vars="VAR_38", variable.name = "egan_id",value.name = "GT")
          
          if(Condition_DEBUG == 1)
          {
            cat("EGAN_FILE_subset.m_1:\n")
            cat(str(EGAN_FILE_subset.m))
            cat("\n")
            cat(str(unique(EGAN_FILE_subset.m$VAR_38)))
            cat("\n")
            cat(str(unique(EGAN_FILE_subset.m$egan_id)))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(EGAN_FILE_subset.m$GT))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(EGAN_FILE_subset.m$GT)))))
            cat("\n")
          }
          
          EGAN_FILE_subset.m$Genotype<-NA
          
          EGAN_FILE_subset.m$Genotype[which(EGAN_FILE_subset.m$GT%in%c('1|0','0|1'))]<-"HET"
          EGAN_FILE_subset.m$Genotype[which(EGAN_FILE_subset.m$GT%in%c('1|1'))]<-"HOM"
          EGAN_FILE_subset.m$Genotype[which(EGAN_FILE_subset.m$GT%in%c('0|0'))]<-"HOM_REF"
          
          EGAN_FILE_subset.m$Genotype<-factor(EGAN_FILE_subset.m$Genotype,
                                              levels=c("HOM_REF","HET","HOM"),
                                              ordered=T)
          
          if(Condition_DEBUG == 1)
          {
            cat("EGAN_FILE_subset.m_2:\n")
            cat(str(EGAN_FILE_subset.m))
            cat("\n")
            cat(str(unique(EGAN_FILE_subset.m$VAR_38)))
            cat("\n")
            cat(str(unique(EGAN_FILE_subset.m$egan_id)))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(EGAN_FILE_subset.m$GT))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(EGAN_FILE_subset.m$GT)))))
            cat("\n")
            cat(sprintf(as.character(names(summary(EGAN_FILE_subset.m$Genotype)))))
            cat("\n")
            cat(sprintf(as.character(summary(EGAN_FILE_subset.m$Genotype))))
            cat("\n")
          }
          
          multiallele_variants_array<-unique(EGAN_FILE_subset.m$VAR_38)
          
          for(iteration_multiallele_variants_array in 1:length(multiallele_variants_array))
          {
            multiallele_variants_array_sel<-multiallele_variants_array[iteration_multiallele_variants_array]
            
            cat("-------------->\t")
            cat(sprintf(as.character(multiallele_variants_array_sel)))
            cat("\n")
            
            EGAN_FILE_subset.m_sel<-EGAN_FILE_subset.m[which(EGAN_FILE_subset.m$VAR_38 == multiallele_variants_array_sel),]
            
            if(Condition_DEBUG == 1)
            {
              cat("EGAN_FILE_subset.m_sel_2:\n")
              cat(str(EGAN_FILE_subset.m_sel))
              cat("\n")
              cat(str(unique(EGAN_FILE_subset.m_sel$VAR_38)))
              cat("\n")
              cat(str(unique(EGAN_FILE_subset.m_sel$egan_id)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(EGAN_FILE_subset.m_sel$GT))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(EGAN_FILE_subset.m_sel$GT)))))
              cat("\n")
              cat(sprintf(as.character(names(summary(EGAN_FILE_subset.m_sel$Genotype)))))
              cat("\n")
              cat(sprintf(as.character(summary(EGAN_FILE_subset.m_sel$Genotype))))
              cat("\n")
            }
            
            df<-merge(EGAN_FILE_subset.m_sel,
                      Thomas_equivalence,
                      by="egan_id")
            if(Condition_DEBUG == 1)
            {
              cat("df:\n")
              cat(str(df))
              cat("\n")
              cat(str(unique(df$VAR_38)))
              cat("\n")
              cat(str(unique(df$egan_id)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(df$GT))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(df$GT)))))
              cat("\n")
              cat(sprintf(as.character(names(summary(df$Genotype)))))
              cat("\n")
              cat(sprintf(as.character(summary(df$Genotype))))
              cat("\n")
            }
            
            df_subset<-unique(df[,c(which(colnames(df) == "VAR_38"),
                                    which(colnames(df) == "sample_id"),
                                    which(colnames(df) == "Genotype"))])
            
            if(Condition_DEBUG == 1)
            {
              cat("df_subset_0:\n")
              cat(str(df_subset))
              cat("\n")
              cat(str(unique(df_subset$VAR_38)))
              cat("\n")
              cat(str(unique(df_subset$sample_id)))
              cat("\n")
              cat(sprintf(as.character(names(summary(df_subset$Genotype)))))
              cat("\n")
              cat(sprintf(as.character(summary(df_subset$Genotype))))
              cat("\n")
            }
            
            df_subset<-merge(df_subset,
                             VAR_and_Proxys,
                             by="VAR_38")
            
            
            if(Condition_DEBUG == 1)
            {
              cat("df_subset_1:\n")
              cat(str(df_subset))
              cat("\n")
              cat(str(unique(df_subset$VAR_38)))
              cat("\n")
              cat(str(unique(df_subset$name)))
              cat("\n")
              cat(str(unique(df_subset$sample_id)))
              cat("\n")
              cat(sprintf(as.character(names(summary(df_subset$Genotype)))))
              cat("\n")
              cat(sprintf(as.character(summary(df_subset$Genotype))))
              cat("\n")
              
              cat("Proxy_file_REMINDER:\n")
              cat(str(Proxy_file))
              cat("\n")
            }
            
            
            if(dim(df_subset)[1] >0)
            {
              
              df_subset_reduced<-unique(df_subset[,c(which(colnames(df_subset) == "sample_id"),
                                                     which(colnames(df_subset) == "Genotype"))])
              
              # if(Condition_DEBUG == 1)
              # {
              cat("df_subset_reduced_0:\n")
              cat(str(df_subset_reduced))
              cat("\n")
              cat(sprintf(as.character(names(summary(df_subset_reduced$Genotype)))))
              cat("\n")
              cat(sprintf(as.character(summary(df_subset_reduced$Genotype))))
              cat("\n")
              
              # }
              
              
              check_RNA_id_df<-unique(Thomas_equivalence[-which(Thomas_equivalence$egan_id%in%EGAN_FILE_subset.m_sel$egan_id),])
              
              # cat("check_RNA_id_df:\n")
              # cat(str(check_RNA_id_df))
              # cat("\n")
              
              
              if(dim(check_RNA_id_df)[1] >0)
              {
                
                
                check_RNA_id_df_subset<-check_RNA_id_df[,c(which(colnames(check_RNA_id_df) == "sample_id"),
                                                           which(colnames(check_RNA_id_df) == "egan_id"))]
                
                # cat("WARNING missing RNA files:\n")
                # cat(str(check_RNA_id_df_subset))
                # cat("\n")
                
                ABSENT_RNA_ids[[i]]<-check_RNA_id_df_subset
                
                
                # setwd(out)
                # 
                # write.table(check_RNA_id_df_subset, file="Absent_EGAN_IDS_in_WGS_vcf.tsv", sep="\t",quote=F, row.names=F)
                # quit(status = 1)
                
              }# dim(check_RNA_id_df)[1] >0
              
              
              
              Proxy_file_VAR<-Proxy_file[which(Proxy_file$VAR%in%df_subset$name),]
              
              
              Proxy_file_Proxy_VAR<-Proxy_file[which(Proxy_file$Proxy_VAR%in%df_subset$name),]
              
              
              
              if(dim(Proxy_file_VAR)[1] >0 & dim(Proxy_file_Proxy_VAR)[1] ==0)
              {
                cat("----Main_Variant--with Proxys---but_not_Proxy_of_other_Main_vars>\n")
                Condition_DEBUG <- 0
                
                # quit(status = 1)
                if(Condition_DEBUG == 1)
                {
                  cat("Proxy_file_VAR_0:\n")
                  cat(str(Proxy_file_VAR))
                  cat("\n")
                  cat(str(unique(Proxy_file_VAR$VAR)))
                  cat("\n")
                  cat(str(unique(Proxy_file_VAR$Proxy_VAR)))
                  cat("\n")
                }
                
                VAR_sel<-unique(Proxy_file_VAR$VAR)
                
                cat("--------------->\t")
                cat(sprintf(as.character(VAR_sel)))
                cat("\n")
                
                path6<-paste(out,VAR_sel,'/', sep='')
                
                # cat("path6\n")
                # cat(sprintf(as.character(path6)))
                # cat("\n")
                
                
                if (file.exists(path6)){
                  
                  
                  
                  
                } else {
                  dir.create(file.path(path6))
                  
                }
                
                ######## covariate specific assignation----
                
                ALL_dB_Thresholded_sel<-ALL_dB_Thresholded[which(ALL_dB_Thresholded$VAR %in% VAR_sel),]
                
                if(Condition_DEBUG == 1)
                {
                  cat("ALL_dB_Thresholded_sel_1\n")
                  cat(str(ALL_dB_Thresholded_sel))
                  cat("\n")
                }
                
                covariate_phenotype_correspondence_sel<-covariate_phenotype_correspondence[which(covariate_phenotype_correspondence$phenotype%in%ALL_dB_Thresholded_sel$phenotype),]
                
                if(Condition_DEBUG == 1)
                {
                  cat("covariate_phenotype_correspondence_sel\n")
                  cat(str(covariate_phenotype_correspondence_sel))
                  cat("\n")
                }
                
                specific_covariates<-unique(covariate_phenotype_correspondence_sel$covariates)
                
                if(Condition_DEBUG == 1)
                {
                  cat("specific_covariates\n")
                  cat(str(specific_covariates))
                  cat("\n")
                }
                
                
                #### Select the variant specific covariates to perform the LM regression ----
                
                
                CELL_covariates<-specific_covariates
                
                CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                      "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                      "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                      "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                      "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")
                
                INTERVAL_covariates_and_PEER_factors_sel<-INTERVAL_covariates_and_PEER_factors[,c(which(colnames(INTERVAL_covariates_and_PEER_factors) == "sample_id"),
                                                                                                  which(colnames(INTERVAL_covariates_and_PEER_factors)%in%c(CONSTITUTIVE_TERMS,CELL_covariates)))]
                
                if(Condition_DEBUG == 1)
                {
                  cat("INTERVAL_covariates_and_PEER_factors_sel\n")
                  cat(str(INTERVAL_covariates_and_PEER_factors_sel))
                  cat("\n")
                }
                
                if(Condition_DEBUG == 1)
                {
                  cat("df_subset_reduced\n")
                  cat(str(df_subset_reduced))
                  cat("\n")
                }
                
                
                INTERVAL_covariates_and_PEER_factors_sel_ADAPTED<-merge(INTERVAL_covariates_and_PEER_factors_sel,
                                                                        df_subset_reduced,
                                                                        by="sample_id")
                
                
                
                if(Condition_DEBUG == 1)
                {
                  cat("INTERVAL_covariates_and_PEER_factors_sel_ADAPTED\n")
                  cat(str(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED))
                  cat("\n")
                }
                
                setwd(path6)
                
                saveRDS(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED, file=paste("INTERVAL_covariates_and_PEER_factors_",VAR_sel,".rds", sep=''))
                
                # quit(status = 1)
                
              }# dim(Proxy_file_VAR)[1] >0 & dim(Proxy_file_Proxy_VAR)[1] ==0
              
              if(dim(Proxy_file_VAR)[1] ==0 & dim(Proxy_file_Proxy_VAR)[1] >0)
              {
                cat("----Proxy_variant ONLY----->\t")
                
                # quit(status = 1)
                
                Condition_DEBUG <- 1
                
                if(Condition_DEBUG == 1)
                {
                  cat("Proxy_file_Proxy_VAR_0:\n")
                  cat(str(Proxy_file_Proxy_VAR))
                  cat("\n")
                  cat(str(unique(Proxy_file_Proxy_VAR$VAR)))
                  cat("\n")
                  cat(str(unique(Proxy_file_Proxy_VAR$Proxy_VAR)))
                  cat("\n")
                }
                
                Proxy_VAR_sel<-unique(Proxy_file_Proxy_VAR$Proxy_VAR)
                
                # cat("--------------->\t")
                cat(sprintf(as.character(Proxy_VAR_sel)))
                cat("\n")
                
                VAR_array<-unique(Proxy_file_Proxy_VAR$VAR)
                
                for(k in 1:length(VAR_array))
                {
                  VAR_array_sel<-VAR_array[k]
                  
                  cat("-----VAR>\t")
                  cat(sprintf(as.character(VAR_array_sel)))
                  cat("\n")
                  
                  # #############################################################################
                  # quit(status = 1)
                  
                  path6<-paste(out,VAR_array_sel,'/', sep='')
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path6)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path6))
                    
                  }
                  
                  path7<-paste(out,VAR_array_sel,'/','Proxys','/', sep='')
                  
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path7)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path7))
                    
                  }
                  
                  path8<-paste(out,VAR_array_sel,'/','Proxys','/',Proxy_VAR_sel,'/', sep='')
                  
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path8)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path8))
                    
                  }
                  
                  ######## covariate specific assignation----
                  
                  ALL_dB_Thresholded_sel<-ALL_dB_Thresholded[which(ALL_dB_Thresholded$VAR %in%Proxy_VAR_sel),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("ALL_dB_Thresholded_sel_1\n")
                    cat(str(ALL_dB_Thresholded_sel))
                    cat("\n")
                  }
                  
                  covariate_phenotype_correspondence_sel<-covariate_phenotype_correspondence[which(covariate_phenotype_correspondence$phenotype%in%ALL_dB_Thresholded_sel$phenotype),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("covariate_phenotype_correspondence_sel\n")
                    cat(str(covariate_phenotype_correspondence_sel))
                    cat("\n")
                  }
                  
                  specific_covariates<-unique(covariate_phenotype_correspondence_sel$covariates)
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("specific_covariates\n")
                    cat(str(specific_covariates))
                    cat("\n")
                  }
                  
                  
                  #### Select the variant specific covariates to perform the LM regression ----
                  
                  
                  CELL_covariates<-specific_covariates
                  
                  CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                        "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                        "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                        "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                        "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")
                  
                  INTERVAL_covariates_and_PEER_factors_sel<-INTERVAL_covariates_and_PEER_factors[,c(which(colnames(INTERVAL_covariates_and_PEER_factors) == "sample_id"),
                                                                                                    which(colnames(INTERVAL_covariates_and_PEER_factors)%in%c(CONSTITUTIVE_TERMS,CELL_covariates)))]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_covariates_and_PEER_factors_sel\n")
                    cat(str(INTERVAL_covariates_and_PEER_factors_sel))
                    cat("\n")
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("df_subset_reduced\n")
                    cat(str(df_subset_reduced))
                    cat("\n")
                  }
                  
                  
                  INTERVAL_covariates_and_PEER_factors_sel_ADAPTED<-merge(INTERVAL_covariates_and_PEER_factors_sel,
                                                                          df_subset_reduced,
                                                                          by="sample_id")
                  
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_covariates_and_PEER_factors_sel_ADAPTED\n")
                    cat(str(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED))
                    cat("\n")
                  }
                  
                  setwd(path8)
                  
                  saveRDS(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED, file=paste("INTERVAL_covariates_and_PEER_factors_",Proxy_VAR_sel,".rds", sep=''))
                  
                  
                  
                }#k in 1:length(VAR_array)
                
              }# dim(Proxy_file_VAR)[1] ==0 & dim(Proxy_file_Proxy_VAR)[1] >0
              
              if(dim(Proxy_file_VAR)[1] >0 & dim(Proxy_file_Proxy_VAR)[1] >0)
              {
                cat("----Main_Variant--with Proxys---AND_Proxy_of_other_Main_vars>\n")
                
                Condition_DEBUG <- 1
                
                if(Condition_DEBUG == 1)
                {
                  cat("Proxy_file_VAR_0:\n")
                  cat(str(Proxy_file_VAR))
                  cat("\n")
                  cat(str(unique(Proxy_file_VAR$VAR)))
                  cat("\n")
                  cat(str(unique(Proxy_file_VAR$Proxy_VAR)))
                  cat("\n")
                }
                
                # ##############################################
                # quit(status = 1)
                
                VAR_sel<-unique(Proxy_file_VAR$VAR)
                
                cat("--------------->\t")
                cat(sprintf(as.character(VAR_sel)))
                cat("\n")
                
                path6<-paste(out,VAR_sel,'/', sep='')
                
                # cat("path6\n")
                # cat(sprintf(as.character(path6)))
                # cat("\n")
                
                
                if (file.exists(path6)){
                  
                  
                  
                  
                } else {
                  dir.create(file.path(path6))
                  
                }
                
                ######## covariate specific assignation----
                
                ALL_dB_Thresholded_sel<-ALL_dB_Thresholded[which(ALL_dB_Thresholded$VAR %in% VAR_sel),]
                
                if(Condition_DEBUG == 1)
                {
                  cat("ALL_dB_Thresholded_sel_1\n")
                  cat(str(ALL_dB_Thresholded_sel))
                  cat("\n")
                }
                
                covariate_phenotype_correspondence_sel<-covariate_phenotype_correspondence[which(covariate_phenotype_correspondence$phenotype%in%ALL_dB_Thresholded_sel$phenotype),]
                
                if(Condition_DEBUG == 1)
                {
                  cat("covariate_phenotype_correspondence_sel\n")
                  cat(str(covariate_phenotype_correspondence_sel))
                  cat("\n")
                }
                
                specific_covariates<-unique(covariate_phenotype_correspondence_sel$covariates)
                
                if(Condition_DEBUG == 1)
                {
                  cat("specific_covariates\n")
                  cat(str(specific_covariates))
                  cat("\n")
                }
                
                
                #### Select the variant specific covariates to perform the LM regression ----
                
                
                CELL_covariates<-specific_covariates
                
                CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                      "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                      "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                      "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                      "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")
                
                INTERVAL_covariates_and_PEER_factors_sel<-INTERVAL_covariates_and_PEER_factors[,c(which(colnames(INTERVAL_covariates_and_PEER_factors) == "sample_id"),
                                                                                                  which(colnames(INTERVAL_covariates_and_PEER_factors)%in%c(CONSTITUTIVE_TERMS,CELL_covariates)))]
                
                if(Condition_DEBUG == 1)
                {
                  cat("INTERVAL_covariates_and_PEER_factors_sel\n")
                  cat(str(INTERVAL_covariates_and_PEER_factors_sel))
                  cat("\n")
                }
                
                if(Condition_DEBUG == 1)
                {
                  cat("df_subset_reduced\n")
                  cat(str(df_subset_reduced))
                  cat("\n")
                }
                
                
                INTERVAL_covariates_and_PEER_factors_sel_ADAPTED<-merge(INTERVAL_covariates_and_PEER_factors_sel,
                                                                        df_subset_reduced,
                                                                        by="sample_id")
                
                
                
                if(Condition_DEBUG == 1)
                {
                  cat("INTERVAL_covariates_and_PEER_factors_sel_ADAPTED\n")
                  cat(str(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED))
                  cat("\n")
                }
                
                setwd(path6)
                
                saveRDS(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED, file=paste("INTERVAL_covariates_and_PEER_factors_",VAR_sel,".rds", sep=''))
                
                # quit(status = 1)
                
                ##########################################################################################################
                
                Proxy_VAR_sel<-unique(Proxy_file_Proxy_VAR$Proxy_VAR)
                
                # cat("--------------->\t")
                cat(sprintf(as.character(Proxy_VAR_sel)))
                cat("\n")
                
                VAR_array<-unique(Proxy_file_Proxy_VAR$VAR)
                
                for(k in 1:length(VAR_array))
                {
                  VAR_array_sel<-VAR_array[k]
                  
                  cat("-----VAR>\t")
                  cat(sprintf(as.character(VAR_array_sel)))
                  cat("\n")
                  
                  # #############################################################################
                  # quit(status = 1)
                  
                  path6<-paste(out,VAR_array_sel,'/', sep='')
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path6)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path6))
                    
                  }
                  
                  path7<-paste(out,VAR_array_sel,'/','Proxys','/', sep='')
                  
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path7)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path7))
                    
                  }
                  
                  path8<-paste(out,VAR_array_sel,'/','Proxys','/',Proxy_VAR_sel,'/', sep='')
                  
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path8)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path8))
                    
                  }
                  
                  ######## covariate specific assignation----
                  
                  ALL_dB_Thresholded_sel<-ALL_dB_Thresholded[which(ALL_dB_Thresholded$VAR %in%Proxy_VAR_sel),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("ALL_dB_Thresholded_sel_1\n")
                    cat(str(ALL_dB_Thresholded_sel))
                    cat("\n")
                  }
                  
                  covariate_phenotype_correspondence_sel<-covariate_phenotype_correspondence[which(covariate_phenotype_correspondence$phenotype%in%ALL_dB_Thresholded_sel$phenotype),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("covariate_phenotype_correspondence_sel\n")
                    cat(str(covariate_phenotype_correspondence_sel))
                    cat("\n")
                  }
                  
                  specific_covariates<-unique(covariate_phenotype_correspondence_sel$covariates)
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("specific_covariates\n")
                    cat(str(specific_covariates))
                    cat("\n")
                  }
                  
                  
                  #### Select the variant specific covariates to perform the LM regression ----
                  
                  
                  CELL_covariates<-specific_covariates
                  
                  CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                        "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                        "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                        "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                        "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")
                  
                  INTERVAL_covariates_and_PEER_factors_sel<-INTERVAL_covariates_and_PEER_factors[,c(which(colnames(INTERVAL_covariates_and_PEER_factors) == "sample_id"),
                                                                                                    which(colnames(INTERVAL_covariates_and_PEER_factors)%in%c(CONSTITUTIVE_TERMS,CELL_covariates)))]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_covariates_and_PEER_factors_sel\n")
                    cat(str(INTERVAL_covariates_and_PEER_factors_sel))
                    cat("\n")
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("df_subset_reduced\n")
                    cat(str(df_subset_reduced))
                    cat("\n")
                  }
                  
                  
                  INTERVAL_covariates_and_PEER_factors_sel_ADAPTED<-merge(INTERVAL_covariates_and_PEER_factors_sel,
                                                                          df_subset_reduced,
                                                                          by="sample_id")
                  
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_covariates_and_PEER_factors_sel_ADAPTED\n")
                    cat(str(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED))
                    cat("\n")
                  }
                  
                  setwd(path8)
                  
                  saveRDS(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED, file=paste("INTERVAL_covariates_and_PEER_factors_",Proxy_VAR_sel,".rds", sep=''))
                  
                  
                  
                }#k in 1:length(VAR_array)
                
              }# dim(Proxy_file_VAR)[1] >0 & dim(Proxy_file_Proxy_VAR)[1] >0
              
              if(dim(Proxy_file_VAR)[1] ==0 & dim(Proxy_file_Proxy_VAR)[1] ==0)
              {
                SELECTED_VARS_sel<-SELECTED_VARS[which(SELECTED_VARS%in%df_subset$name)]
                
                cat("SELECTED_VARS_sel\n")
                cat(sprintf(as.character(SELECTED_VARS_sel)))
                cat("\n")
                
                if(length(SELECTED_VARS_sel) >0)
                {
                  cat("----Main_Variant--without_Proxys--->\n")
                  Condition_DEBUG <- 0
                  
                  # quit(status = 1)
                  
                  path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
                  
                  # cat("path6\n")
                  # cat(sprintf(as.character(path6)))
                  # cat("\n")
                  
                  
                  if (file.exists(path6)){
                    
                    
                    
                    
                  } else {
                    dir.create(file.path(path6))
                    
                  }
                  
                  ######## covariate specific assignation----
                  
                  ALL_dB_Thresholded_sel<-ALL_dB_Thresholded[which(ALL_dB_Thresholded$VAR %in% SELECTED_VARS_sel),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("ALL_dB_Thresholded_sel_1\n")
                    cat(str(ALL_dB_Thresholded_sel))
                    cat("\n")
                  }
                  
                  covariate_phenotype_correspondence_sel<-covariate_phenotype_correspondence[which(covariate_phenotype_correspondence$phenotype%in%ALL_dB_Thresholded_sel$phenotype),]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("covariate_phenotype_correspondence_sel\n")
                    cat(str(covariate_phenotype_correspondence_sel))
                    cat("\n")
                  }
                  
                  specific_covariates<-unique(covariate_phenotype_correspondence_sel$covariates)
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("specific_covariates\n")
                    cat(str(specific_covariates))
                    cat("\n")
                  }
                  
                  
                  #### Select the variant specific covariates to perform the LM regression ----
                  
                  
                  CELL_covariates<-specific_covariates
                  
                  CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                        "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                        "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                        "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                        "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")
                  
                  INTERVAL_covariates_and_PEER_factors_sel<-INTERVAL_covariates_and_PEER_factors[,c(which(colnames(INTERVAL_covariates_and_PEER_factors) == "sample_id"),
                                                                                                    which(colnames(INTERVAL_covariates_and_PEER_factors)%in%c(CONSTITUTIVE_TERMS,CELL_covariates)))]
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_covariates_and_PEER_factors_sel\n")
                    cat(str(INTERVAL_covariates_and_PEER_factors_sel))
                    cat("\n")
                  }
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("df_subset_reduced\n")
                    cat(str(df_subset_reduced))
                    cat("\n")
                  }
                  
                  
                  INTERVAL_covariates_and_PEER_factors_sel_ADAPTED<-merge(INTERVAL_covariates_and_PEER_factors_sel,
                                                                          df_subset_reduced,
                                                                          by="sample_id")
                  
                  
                  
                  if(Condition_DEBUG == 1)
                  {
                    cat("INTERVAL_covariates_and_PEER_factors_sel_ADAPTED\n")
                    cat(str(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED))
                    cat("\n")
                  }
                  
                  setwd(path6)
                  
                  saveRDS(INTERVAL_covariates_and_PEER_factors_sel_ADAPTED, file=paste("INTERVAL_covariates_and_PEER_factors_",SELECTED_VARS_sel,".rds", sep=''))
                  
                  
                }else{
                  
                  cat("WARNING_NO_CORRESPONDENCE>\n")
                  
                  quit(status = 1)
                  
                }# length(SELECTED_VARS_sel) >0
                
              }#dim(Proxy_file_VAR)[1] ==0 & dim(Proxy_file_Proxy_VAR)[1] ==0
            }# dim(df_subset)[1] >0
          }#iteration_multiallele_variants_array in 1:length(multiallele_variants_array)
          
          # check_egan_id<-EGAN_FILE_subset.m$egan_id[which(EGAN_FILE_subset.m$egan_id%in%Thomas_equivalence$egan_id)]
          
         
        }#dim(EGAN_FILE)[1] >0
        
      }#n_lines > indx_header
        
        
    }#file.exists(file_sel)
    
  }#i file_list_sel
  
  
  ### Print the list of missing RNA samples with EGAN IDs from WES
  
}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----<

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
    make_option(c("--ALL_dB"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--covariate_phenotype_correspondence"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--INTERVAL_covariates_and_PEER_factors"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Thomas_equivalence"), type="character", default=NULL, 
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
  
  Data_wrangling(opt)
  
}


###########################################################################

system.time( main() )
