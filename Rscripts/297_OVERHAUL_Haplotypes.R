

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
  
  #### Selected_vars ----
  
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS))
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
  
  
  
  SELECTED_VARS_UPDATED = SELECTED_VARS[-which(SELECTED_VARS%in%ABSENT_WGS_RNA$VAR)]
  
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
  
  # SELECTED_VARS_UPDATED<-"chr2_144162105_A_G"
  
  #### LOOP ----
  
  Condition_DEBUG <- 0
  
  for(i in 1:length(SELECTED_VARS_UPDATED))
  {
    Condition_DEBUG <- 0
    
    SELECTED_VARS_UPDATED_sel<-SELECTED_VARS_UPDATED[i]
    
    
    #chr9_135920196_C_T
    #chr8_130641322_C_T
    
    if(SELECTED_VARS_UPDATED_sel == "chr2_144162105_A_G")
    {
      Condition_DEBUG <- 1
    }else{
      
      Condition_DEBUG <- 0
    }
    
    
    
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
        
        path6<-paste(out,SELECTED_VARS_UPDATED_sel,'/', sep='')
        
        # cat("path6\n")
        # cat(sprintf(as.character(path6)))
        # cat("\n")
        
        setwd(path6)
        
        filename_rds_VAR<-paste("INTERVAL_covariates_and_PEER_factors_",SELECTED_VARS_UPDATED_sel,".rds", sep='')
          
        
        
        if (file.exists(filename_rds_VAR)){
         
          cat("Hello_world_1\n")
          
          PEER_G_VAR<-readRDS(file=filename_rds_VAR)
          
          if(Condition_DEBUG == 1)
          {
            cat("PEER_G_VAR\n")
            cat(str(PEER_G_VAR))
            cat("\n")
          }
          
          
        } else {
          
          cat("WARNING_ABSENT_PEER_G_MAIN_VAR\n")
          quit(status = 1)
        }
        
        path8<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Proxys','/',Proxy_array_sel,'/', sep='')
        
        if (file.exists(path8)){
          
          setwd(path8)
          
          filename_rds_Proxy_VAR<-paste("INTERVAL_covariates_and_PEER_factors_",Proxy_array_sel,".rds", sep='')
          
          
          
          if (file.exists(filename_rds_Proxy_VAR)){
            
            PEER_G_Proxy_VAR<-readRDS(file=filename_rds_Proxy_VAR)
            
            if(Condition_DEBUG == 1)
            {
              cat("PEER_G_Proxy_VAR\n")
              cat(str(PEER_G_Proxy_VAR))
              cat("\n")
            }
            
            
          } else {
            
            cat("WARNING_ABSENT_PEER_G_PROXY_VAR\n")
            quit(status = 1)
          }
          
          
          ## Select samples in common and covariates to the maximum
          
          CONSTITUTIVE_TERMS<-c("age_RNA","BMI","Conc_ng_ul","RIN","RawReadDepth","Season_Winter","Season_Autumn","Season_Spring","Season_Summer","sequencingBatch_1","sequencingBatch_2","sequencingBatch_3","sequencingBatch_4","sequencingBatch_5",
                                "sequencingBatch_6","sequencingBatch_7","sequencingBatch_8","sequencingBatch_9","sequencingBatch_10","sequencingBatch_11","sequencingBatch_12","sequencingBatch_13","sequencingBatch_14","sequencingBatch_15",
                                "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                                "PEER1","PEER2","PEER3","PEER4","PEER5","PEER6","PEER7","PEER8","PEER9","PEER10","PEER11","PEER12","PEER13","PEER14","PEER15","PEER16","PEER17","PEER18","PEER19","PEER20","PEER21","PEER22","PEER23","PEER24",
                                "PEER25","PEER26","PEER27","PEER28","PEER29","PEER30","PEER31","PEER32","PEER33","PEER34","PEER35","sex")
          
          genotype<-"Genotype"
          sample_id<-"sample_id"
          
          cell_specific_covariates_VAR<-colnames(PEER_G_VAR)[-which(colnames(PEER_G_VAR)%in%c(CONSTITUTIVE_TERMS,genotype,sample_id))]
          
          if(Condition_DEBUG == 1)
          {
            cat("cell_specific_covariates_VAR\n")
            cat(str(cell_specific_covariates_VAR))
            cat("\n")
          }
          
          cell_specific_covariates_Proxy_VAR<-colnames(PEER_G_Proxy_VAR)[-which(colnames(PEER_G_Proxy_VAR)%in%c(CONSTITUTIVE_TERMS,genotype,sample_id))]
          
          if(Condition_DEBUG == 1)
          {
            cat("cell_specific_covariates_Proxy_VAR\n")
            cat(str(cell_specific_covariates_Proxy_VAR))
            cat("\n")
          }
          
          Haplotype_specific_covariates<-unique(c(cell_specific_covariates_VAR,cell_specific_covariates_Proxy_VAR))
          
          if(Condition_DEBUG == 1)
          {
            cat("Haplotype_specific_covariates\n")
            cat(str(Haplotype_specific_covariates))
            cat("\n")
          }
          
          
          Haplotype_colnames<-unique(c(genotype,sample_id,CONSTITUTIVE_TERMS,Haplotype_specific_covariates))
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Haplotype_colnames\n")
            cat(str(Haplotype_colnames))
            cat("\n")
          }
          
          
          intersect_samples<-PEER_G_VAR$sample_id[which(PEER_G_VAR$sample_id%in%PEER_G_Proxy_VAR$sample_id)]
          
          if(Condition_DEBUG == 1)
          {
            cat("intersect_samples\n")
            cat(str(intersect_samples))
            cat("\n")
          }
          
          PEER_G_VAR_Haplotype<-PEER_G_VAR[,which(colnames(PEER_G_VAR)%in%Haplotype_colnames)]
          
          
          
          # PEER_G_VAR_Haplotype<-PEER_G_VAR_Haplotype[,Haplotype_colnames]
          
          if(Condition_DEBUG == 1)
          {
            # cat("PEER_G_VAR_Haplotype_0\n")
            # cat(str(PEER_G_VAR_Haplotype))
            # cat("\n")
          }
          
          PEER_G_VAR_Haplotype<-PEER_G_VAR_Haplotype[which(PEER_G_VAR_Haplotype$sample_id%in%intersect_samples),]
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("PEER_G_VAR_Haplotype_1\n")
          #   cat(str(PEER_G_VAR_Haplotype))
          #   cat("\n")
          # }
          
          PEER_G_Proxy_VAR_Haplotype<-PEER_G_Proxy_VAR[,which(colnames(PEER_G_Proxy_VAR)%in%Haplotype_colnames)]
          
          
          
          # PEER_G_Proxy_VAR_Haplotype<-PEER_G_Proxy_VAR_Haplotype[,Haplotype_colnames]
          
          if(Condition_DEBUG == 1)
          {
            # cat("PEER_G_Proxy_VAR_Haplotype_0\n")
            # cat(str(PEER_G_Proxy_VAR_Haplotype))
            # cat("\n")
          }
          
          PEER_G_Proxy_VAR_Haplotype<-PEER_G_Proxy_VAR_Haplotype[which(PEER_G_Proxy_VAR_Haplotype$sample_id%in%intersect_samples),]
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("PEER_G_Proxy_VAR_Haplotype_1\n")
          #   cat(str(PEER_G_Proxy_VAR_Haplotype))
          #   cat("\n")
          # }
          
          
          common_colnames_minus_Genotype<-colnames(PEER_G_VAR_Haplotype)[which(colnames(PEER_G_VAR_Haplotype)%in%colnames(PEER_G_Proxy_VAR_Haplotype))]
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("common_colnames_minus_Genotype_0\n")
          #   cat(str(common_colnames_minus_Genotype))
          #   cat("\n")
          # }
          
          common_colnames_minus_Genotype<-common_colnames_minus_Genotype[-which(common_colnames_minus_Genotype%in%genotype)]
          
          if(Condition_DEBUG == 1)
          {
            cat("common_colnames_minus_Genotype_1\n")
            cat(str(common_colnames_minus_Genotype))
            cat("\n")
          }
          
          Haplotype_PEER_G<-merge(PEER_G_VAR_Haplotype,
                                  PEER_G_Proxy_VAR_Haplotype,
                                  by=common_colnames_minus_Genotype)
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("Haplotype_PEER_G_0\n")
          #   cat(str(Haplotype_PEER_G))
          #   cat("\n")
          # }
          
          colnames(Haplotype_PEER_G)[which(colnames(Haplotype_PEER_G) == "Genotype.x")]<-"Genotype_VAR"
          colnames(Haplotype_PEER_G)[which(colnames(Haplotype_PEER_G) == "Genotype.y")]<-"Genotype_Proxy_VAR"
          
          indx.first<-unique(c(which(colnames(Haplotype_PEER_G) == "Genotype_VAR"),
                               which(colnames(Haplotype_PEER_G) == "Genotype_Proxy_VAR")))
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("indx.first_0\n")
          #   cat(str(indx.first))
          #   cat("\n")
          # }
          
          indx.second<-which(colnames(Haplotype_PEER_G) == "sample_id")
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("indx.second_0\n")
          #   cat(str(indx.second))
          #   cat("\n")
          # }
          
          indx.REST<-which(colnames(Haplotype_PEER_G) != "Genotype_VAR" &
                             colnames(Haplotype_PEER_G) != "Genotype_Proxy_VAR" &
                             colnames(Haplotype_PEER_G) != "sample_id")
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("indx.REST_0\n")
          #   cat(str(indx.REST))
          #   cat("\n")
          # }
          
          order_DEF<-c(indx.first,indx.second,indx.REST)
          
          if(Condition_DEBUG == 1)
          {
            cat("order_DEF_1\n")
            cat(str(order_DEF))
            cat("\n")
          }
          
          Haplotype_PEER_G<-Haplotype_PEER_G[,order_DEF]
          
          if(Condition_DEBUG == 1)
          {
            cat("Haplotype_PEER_G_1\n")
            cat(str(Haplotype_PEER_G))
            cat("\n")
          }
          
          Haplotype_PEER_G$Haplotype_interaction<-interaction(Haplotype_PEER_G$Genotype_VAR,Haplotype_PEER_G$Genotype_Proxy_VAR, sep="|",lex.order = T)
          
          # Haplotype_PEER_G$Haplotype_interaction<-droplevels(Haplotype_PEER_G$Haplotype_interaction)
          
          if(Condition_DEBUG == 1)
          {
            # cat("Haplotype_PEER_G_2\n")
            # cat(str(Haplotype_PEER_G))
            # cat("\n")
            # cat(sprintf(as.character(names(summary(Haplotype_PEER_G$Haplotype_interaction)))))
            # cat("\n")
            # cat(sprintf(as.character(summary(Haplotype_PEER_G$Haplotype_interaction))))
            # cat("\n")
          }
          
          Haplotype_PEER_G$Haplotype<-"NA"
          
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HOM_REF|HOM_REF')]<-'HOM_REF'
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HET|HOM_REF')]<-'HET|HOM_REF'
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HOM_REF|HET')]<-'HOM_REF|HET'
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HET|HET')]<-'HET|HET'
          
          
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HOM|HOM_REF')]<-'HOM|HOM_REF'
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HOM_REF|HOM')]<-'HOM_REF|HOM'
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HOM|HOM')]<-'HOM|HOM'
          
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HOM|HET')]<-'HOM|HET'
          Haplotype_PEER_G$Haplotype[which(Haplotype_PEER_G$Haplotype_interaction == 'HET|HOM')]<-'HET|HOM'
          
          Haplotype_PEER_G$Haplotype<-factor(Haplotype_PEER_G$Haplotype,
                                             levels=c('HOM_REF','HET|HOM_REF','HOM_REF|HET','HET|HET','HOM|HOM_REF','HOM_REF|HOM','HOM|HET','HET|HOM','HOM|HOM'),
                                             ordered=T)
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Haplotype_PEER_G_3\n")
            cat(str(Haplotype_PEER_G))
            cat("\n")
            cat(sprintf(as.character(names(summary(Haplotype_PEER_G$Genotype_VAR)))))
            cat("\n")
            cat(sprintf(as.character(summary(Haplotype_PEER_G$Genotype_VAR))))
            cat("\n")
            cat(sprintf(as.character(names(summary(Haplotype_PEER_G$Genotype_Proxy_VAR)))))
            cat("\n")
            cat(sprintf(as.character(summary(Haplotype_PEER_G$Genotype_Proxy_VAR))))
            cat("\n")
            
            cat(sprintf(as.character(names(summary(Haplotype_PEER_G$Haplotype_interaction)))))
            cat("\n")
            cat(sprintf(as.character(summary(Haplotype_PEER_G$Haplotype_interaction))))
            cat("\n")
            cat(sprintf(as.character(names(summary(Haplotype_PEER_G$Haplotype)))))
            cat("\n")
            cat(sprintf(as.character(summary(Haplotype_PEER_G$Haplotype))))
            cat("\n")
          }
          
          indx.dep<-c(which(colnames(Haplotype_PEER_G) == "Genotype_VAR"),which(colnames(Haplotype_PEER_G) == "Genotype_Proxy_VAR"),which(colnames(Haplotype_PEER_G) == "Haplotype_interaction"))
          
          Haplotype_PEER_G<-Haplotype_PEER_G[,-indx.dep]
          
          indx.reorder<-c(which(colnames(Haplotype_PEER_G) == "Haplotype"),which(colnames(Haplotype_PEER_G) != "Haplotype"))
          
          Haplotype_PEER_G<-Haplotype_PEER_G[,indx.reorder]
          
          if(Condition_DEBUG == 1)
          {
            cat("Haplotype_PEER_G_4\n")
            cat(str(Haplotype_PEER_G))
            cat("\n")
            cat(sprintf(as.character(names(summary(Haplotype_PEER_G$Haplotype)))))
            cat("\n")
            cat(sprintf(as.character(summary(Haplotype_PEER_G$Haplotype))))
            cat("\n")
          }
          
          
          
          path9<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/', sep='')
          
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path9)){
            
            
            
            
          } else {
            dir.create(file.path(path9))
            
          }
          
          path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
          
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path10)){
            
            
            
            
          } else {
            dir.create(file.path(path10))
            
          }
          
          setwd(path10)
          
          saveRDS(Haplotype_PEER_G, file=paste("INTERVAL_covariates_and_PEER_factors_Haplotype_",paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),".rds", sep=''))
          
          
          
          # ############################################################################################################
          # quit(status = 1)
          
          
        }else{
          
          cat("WARNING_ABSENT_PEER_G_Proxy_VAR\n")
          quit(status = 1)
          
        }#file.exists(path8)
      }#k in 1:length(Proxy_array)
      # ############################################################################################################
      # quit(status = 1)
      
    }#dim(Proxy_file_UPDATED_sel)[1] >0
  }#i in 1:length(SELECTED_VARS_UPDATED)
  
  # ############################################################################################################
  # quit(status = 1)

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
