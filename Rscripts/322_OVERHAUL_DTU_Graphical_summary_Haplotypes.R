

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
  
  
  
  SELECTED_VARS  = unlist(strsplit(opt$SELECTED_VARS , split=","))
  
  cat("SELECTED_VARS_\n")
  cat(str(SELECTED_VARS ))
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
  
  indx.UPDATE<-which(SELECTED_VARS %in%ABSENT_WGS_RNA$VAR)
  
  if(length(indx.UPDATE) >0)
  {
    SELECTED_VARS_UPDATED = SELECTED_VARS [-indx.UPDATE]
    
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
  
  
  ##### LOOP TO READ ALL VARIABLES -----
  
  
  Condition_DEBUG <- 0
  
  # SELECTED_VARS_UPDATED<-"chr9_135920196_C_T"
  
  
  if(length(SELECTED_VARS_UPDATED) >0)
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
    #### Read transposed expression ----
    
    
    
    Transposed_Isoform_Expression<-readRDS(file=opt$Transposed_Isoform_Expression)
    colnames(Transposed_Isoform_Expression)[which(colnames(Transposed_Isoform_Expression) == "Sample_id")]<-"sample_id"
    
    
   
    # quit(status=1)
    
    
    ### Read Results_LogLM----
    
    
    
    Results_LogLM<-as.data.frame(fread(file=opt$Results_LogLM, sep="\t", header=T) , stringsAsFactors=T)
    
    
    cat("Results_LogLM\n")
    cat(str(Results_LogLM))
    cat("\n")
    cat(str(unique(Results_LogLM$VAR)))
    cat("\n")
    
   
    # quit(status = 1)
    
    
    ### size scale ----
    
    minuslospval_vector_INTERVAL<-unique(c(Results_LogLM$CIS_gene_minuslogpvalue,Results_LogLM$Block_PCHiC_minuslogpvalue))
    
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
    
    ##### MAIN LOOP------
    
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
        
        
        for(k in 1:length(Proxy_array))
        {
          Proxy_array_sel<-Proxy_array[k]
          
          
          cat("Proxy_VAR---------->\t")
          cat(sprintf(as.character(Proxy_array_sel)))
          cat("\n")
          
          # if(SELECTED_VARS_UPDATED_sel =="chr9_135920196_C_T" & Proxy_array_sel == "chr9_135864513_C_T")
          # if(SELECTED_VARS_UPDATED_sel =="chr16_67250992_C_T" & Proxy_array_sel == "chr16_66895531_G_A")
          # {
          #   Condition_DEBUG <- 1
          # }else{
          #   Condition_DEBUG <- 0
          # }
          
          
          
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
              
              # if(Condition_DEBUG == 1)
              # {
              #   cat("Haplotype_PEER_G\n")
              #   cat(str(Haplotype_PEER_G))
              #   cat("\n")
              # }
              
              Haplotype_PEER_G_HET_haplotypes<-droplevels(Haplotype_PEER_G[(as.numeric(Haplotype_PEER_G$Haplotype) <= 4),])
              
              Haplotype_PEER_G_HET_haplotypes$Haplotype<-factor(Haplotype_PEER_G_HET_haplotypes$Haplotype,
                                                                levels=c('HOM_REF','HET|HOM_REF','HOM_REF|HET','HET|HET'),
                                                                ordered=T)
              
              # if(Condition_DEBUG == 1)
              # {
              #   cat("Haplotype_PEER_G_HET_haplotypes\n")
              #   cat(str(Haplotype_PEER_G_HET_haplotypes))
              #   cat("\n")
              # }
              
              
              
              Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF<-droplevels(Haplotype_PEER_G_HET_haplotypes[(as.numeric(Haplotype_PEER_G_HET_haplotypes$Haplotype) > 1),])
              
              
              # if(Condition_DEBUG == 1)
              # {
              #   cat("Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF\n")
              #   cat(str(Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF))
              #   cat("\n")
              # }
              
              samples_NO_HOM_REF<-Haplotype_PEER_G_HET_haplotypes_NO_HOM_REF$sample_id
              
              if(Condition_DEBUG == 1)
              {
                cat("samples_NO_HOM_REF\n")
                cat(str(samples_NO_HOM_REF))
                cat("\n")
              }
              
              Haplotype_PEER_G_HET_haplotypes_HOM_REF<-droplevels(Haplotype_PEER_G_HET_haplotypes[(as.numeric(Haplotype_PEER_G_HET_haplotypes$Haplotype) == 1),])
              
              # if(Condition_DEBUG == 1)
              # {
              #   cat("Haplotype_PEER_G_HET_haplotypes_HOM_REF\n")
              #   cat(str(Haplotype_PEER_G_HET_haplotypes_HOM_REF))
              #   cat("\n")
              # }
              
              samples_HOM_REF<-Haplotype_PEER_G_HET_haplotypes_HOM_REF$sample_id
              
              if(Condition_DEBUG == 1)
              {
                cat("samples_HOM_REF\n")
                cat(str(samples_HOM_REF))
                cat("\n")
              }
              
              Results_LogLM_sel<-Results_LogLM[which(Results_LogLM$VAR == SELECTED_VARS_UPDATED_sel &
                                                                   Results_LogLM$Proxy_VAR == Proxy_array_sel),]
              
              if(Condition_DEBUG == 1)
              {
                cat("Results_LogLM_sel\n")
                cat(str(Results_LogLM_sel))
                cat("\n")
                cat(str(unique(Results_LogLM_sel$VAR)))
                cat("\n")
                cat(str(unique(Results_LogLM_sel$Proxy_VAR)))
                cat("\n")
              }
              
              ENSG_array<-unique(c(Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$CIS_gene_minuslogpvalue)],
                                   Results_LogLM_sel$ensembl_gene_id[!is.na(Results_LogLM_sel$Block_PCHiC_minuslogpvalue)]))
              
              
              # ENSG_array<-"ENSG00000140939"
              ####################################################################################################################
              
              if(length(ENSG_array) >0)
              {
                
                for(z in 1:length(ENSG_array))
                {
                  
                  ENSG_array_sel<-ENSG_array[z]
                  
                  cat("---------------->\t")
                  cat(sprintf(as.character(ENSG_array_sel)))
                  cat("\t")
                  
                  Results_LogLM_sel_ENSG_sel<-Results_LogLM_sel[which(Results_LogLM_sel$ensembl_gene_id%in%ENSG_array_sel),]
                  
                  
                  if(dim(Results_LogLM_sel_ENSG_sel)[1] >0)
                  {
                    if(Condition_DEBUG == 1)
                    {
                      cat("Results_LogLM_sel_ENSG_sel_\n")
                      cat(str(Results_LogLM_sel_ENSG_sel))
                      cat("\n")
                      cat(str(unique(Results_LogLM_sel_ENSG_sel$ensembl_gene_id)))
                      cat("\n")
                    }
                    
                    HGNC_sel<-unique(Results_LogLM_sel_ENSG_sel$HGNC)
                    ENST_array<-unique(Results_LogLM_sel_ENSG_sel$transcript_id)
                    
                    
                    if(length(ENST_array) >0)
                    {
                      if(Condition_DEBUG == 1)
                      {
                        cat("ENST_array_\n")
                        cat(str(ENST_array))
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
                      
                      DTU_gene_df<-band.genes<-data.frame(matrix(ncol=7,nrow=dim(Results_LogLM_sel_ENSG_sel)[1], 
                                                                 dimnames=list(NULL, c("VAR", "Proxy_VAR","comparison",
                                                                                       "ensembl_gene_id","transcript_id",
                                                                                       "Significance",
                                                                                       "RNASeq_source"))),
                                                          stringsAsFactors = F)
                      
                      DTU_gene_df$VAR<-Results_LogLM_sel_ENSG_sel$VAR
                      DTU_gene_df$Proxy_VAR<-Results_LogLM_sel_ENSG_sel$Proxy_VAR
                      DTU_gene_df$comparison<-Results_LogLM_sel_ENSG_sel$comparison
                      DTU_gene_df$ensembl_gene_id<-Results_LogLM_sel_ENSG_sel$ensembl_gene_id
                      DTU_gene_df$transcript_id<-Results_LogLM_sel_ENSG_sel$transcript_id
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("DTU_gene_df_0\n")
                        cat(str(DTU_gene_df))
                        cat("\n")
                        cat(str(unique(DTU_gene_df$ensembl_gene_id)))
                        cat("\n")
                      }
                      
                      
                      
                      if(dim(CIS_gene_INTERVAL)[1] >0)
                      {
                        Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "comparison"),
                                                                             which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                             which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                             which(colnames(Block_genes_INTERVAL) == "transcript_id"),
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
                                                                       which(colnames(CIS_gene_INTERVAL) == "transcript_id"),
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
                        
                        DTU_gene_df<-merge(DTU_gene_df,
                                           Bind_CIS_Block_INTERVAL,
                                           by=c("comparison","ensembl_gene_id","transcript_id"),
                                           all=T)
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("DTU_gene_df_1\n")
                          cat(str(DTU_gene_df))
                          cat("\n")
                          cat(str(unique(DTU_gene_df$ensembl_gene_id)))
                          cat("\n")
                          cat(str(unique(DTU_gene_df$transcript_id)))
                          cat("\n")
                        }
                        
                        # quit(status = 1)
                        
                      }else{
                        
                        Block_genes_INTERVAL_subset<-Block_genes_INTERVAL[,c(which(colnames(Block_genes_INTERVAL) == "comparison"),
                                                                             which(colnames(Block_genes_INTERVAL) == "ensembl_gene_id"),
                                                                             which(colnames(Block_genes_INTERVAL) == "HGNC"),
                                                                             which(colnames(Block_genes_INTERVAL) == "transcript_id"),
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
                        
                        DTU_gene_df<-merge(DTU_gene_df,
                                           Block_genes_INTERVAL_subset,
                                           by=c("comparison","ensembl_gene_id","transcript_id"),
                                           all=T)
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("DTU_gene_df_1\n")
                          cat(str(DTU_gene_df))
                          cat("\n")
                          cat(str(unique(DTU_gene_df$ensembl_gene_id)))
                          cat("\n")
                          cat(str(unique(DTU_gene_df$transcript_id)))
                          cat("\n")
                        }
                        
                      }# dim(CIS_gene_INTERVAL)[1] >0
                      
                      
                      ### Dot plot Graph ----
                      
                      
                      
                      levels_transcript_id<-unique(DTU_gene_df$transcript_id)
                      
                      DTU_gene_df$transcript_id<-factor(DTU_gene_df$transcript_id,
                                                        levels=levels_transcript_id,
                                                        ordered = T)
                      
                      #### Significance
                      
                      DTU_gene_df$Significance[which(DTU_gene_df$ajusted.minuslogpvalue_Haplotypes >= 1.3)]<-"YES"
                      DTU_gene_df$Significance[which(DTU_gene_df$ajusted.minuslogpvalue_Haplotypes < 1.3)]<-"NO"
                      
                      DTU_gene_df$Significance<-factor(DTU_gene_df$Significance,
                                                       levels=c("NO","YES"),
                                                       ordered=T)
                      if(Condition_DEBUG == 1)
                      {
                        cat("DTU_gene_df$Significance\n")
                        cat(sprintf(as.character(names(summary(DTU_gene_df$Significance)))))
                        cat("\n")
                        cat(sprintf(as.character(summary(DTU_gene_df$Significance))))
                        cat("\n")
                      }
                      
                      #### RNASeq_source
                      
                      DTU_gene_df$RNASeq_source<-"Whole blood"
                      
                      
                      DTU_gene_df$RNASeq_source<-factor(DTU_gene_df$RNASeq_source,
                                                        levels=c("Whole blood"),
                                                        ordered=T)
                      if(Condition_DEBUG == 1)
                      {
                        cat("DTU_gene_df$RNASeq_source\n")
                        cat(sprintf(as.character(names(summary(DTU_gene_df$RNASeq_source)))))
                        cat("\n")
                        cat(sprintf(as.character(summary(DTU_gene_df$RNASeq_source))))
                        cat("\n")
                      }
                      
                      
                      ### update break.size
                      
                      local_minuslogpval_max<-max(DTU_gene_df$ajusted.minuslogpvalue_Haplotypes[!is.na(DTU_gene_df$ajusted.minuslogpvalue_Haplotypes)])
                      
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
                      
                      
                      
                      
                      
                      ### Use Beta from the full model of logRatio
                      
                      ### Labels beta
                      
                      indx.finite<-is.finite(DTU_gene_df$coefficient_Haplotypes)
                      
                      check.finite<-sum(indx.finite)
                      
                      if(check.finite >0)
                      {
                        A<-summary(DTU_gene_df$coefficient_Haplotypes[indx.finite])
                        
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
                      
                      ### Tappas
                      
                      DTU_gene_df<-merge(DTU_gene_df,
                                         Tappas_gff_transcript,
                                         by="transcript_id",
                                         all.x=T)
                      
                      DTU_gene_df$transcript_id<-factor(DTU_gene_df$transcript_id,
                                                        levels=rev(levels_transcript_id),
                                                        ordered=T)
                      
                      DTU_gene_df$Significance<-"NO"
                      DTU_gene_df$Significance[which(DTU_gene_df$ajusted.minuslogpvalue_Haplotypes >= 1.3)]<-"YES"
                      
                      DTU_gene_df$Significance<-factor(DTU_gene_df$Significance,
                                                       levels=c("NO","YES"),
                                                       ordered=T)
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("DTU_gene_df_for_dot_plot_1:\t")
                        cat(str(DTU_gene_df))
                        cat("\n")
                        
                      }
                      
                      DTU_gene_df<-merge(DTU_gene_df,
                                         Tappas_gff_CDS,
                                         by="transcript_id",
                                         all.x=T)
                      
                      Tappas_gff_miRNA_Binding_sel<-Tappas_gff_miRNA_Binding[which(Tappas_gff_miRNA_Binding$transcript_id%in%DTU_gene_df$transcript_id),]
                      
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
                      
                      Tappas_gff_DOMAIN_sel<-Tappas_gff_DOMAIN[which(Tappas_gff_DOMAIN$transcript_id%in%DTU_gene_df$transcript_id),]
                      
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
                      
                      Tappas_gff_PTM_sel<-Tappas_gff_PTM[which(Tappas_gff_PTM$transcript_id%in%DTU_gene_df$transcript_id),]
                      
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
                      
                      
                      ### Dot plot Graph
                      
                      DTU_gene_df<-DTU_gene_df[order(DTU_gene_df$transcript_id),]
                      
                      DTU_gene_df$comparison<-factor(DTU_gene_df$comparison,
                                                     levels=c('HOM_REF__HET|HOM_REF',
                                                              'HOM_REF__HOM_REF|HET',
                                                              'HOM_REF__HET|HET'),
                                                     ordered=T)
                      
                      if(Condition_DEBUG == 1)
                      {
                        cat("DTU_gene_df_for_dot_plot_2:\t")
                        cat(str(DTU_gene_df))
                        cat("\n")
                        cat(sprintf(as.character(names(summary(DTU_gene_df$comparison)))))
                        cat("\n")
                        cat(sprintf(as.character(summary(DTU_gene_df$comparison))))
                        cat("\n")
                      }
                      
                      dotplot_INTERVAL_DTU<-DTU_gene_df %>%
                        mutate(myaxis = paste0(transcript_id, "\n", BIOTYPE,"\n",UniProt_ID)) %>%
                        mutate(myaxis=fct_reorder(myaxis,as.numeric(DTU_gene_df$transcript_id))) %>%
                        ggplot(aes(y=myaxis,
                                   x=comparison)) +
                        geom_point(aes(color=Significance,fill=coefficient_Haplotypes,size=ajusted.minuslogpvalue_Haplotypes), stroke=1, shape=21)+
                        scale_color_manual(values=c("gray","black"),name='Significant', drop=F)+
                        scale_size(range = c(0,20), name='-log10pval',
                                   breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)])) +
                        scale_fill_gradient2(
                          low = "blue", 
                          mid = "white", 
                          high = "red", 
                          midpoint = 0,
                          breaks=breaks.Beta_INTERVAL,labels=labels.Beta_INTERVAL,
                          limits=c(breaks.Beta_INTERVAL[1]-0.01,breaks.Beta_INTERVAL[length(breaks.Beta_INTERVAL)]+0.01),name=paste('Effect size','Z-score',sep="\n"),na.value = "gray")+
                        scale_y_discrete(name=NULL, drop=F)+
                        scale_x_discrete(name=NULL, drop=F)+
                        theme_classic()+
                        ggeasy::easy_center_title()
                      
                      
                      dotplot_INTERVAL_DTU<-dotplot_INTERVAL_DTU+
                        facet_grid(cols = vars(DTU_gene_df$comparison), scales='free_x', space='free_x') +
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
                      
                      
                      
                      
                      graph_DEF<-plot_grid(dotplot_INTERVAL_DTU,
                                           nrow = 1,
                                           ncol = 1)
                      
                      
                      
                      
                      svgname<-paste("Dotplot_DTU_",HGNC_sel,".svg",sep='')
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
                      
                      #### Violin plot residuals ----
                      
                      ### read residuals files
                      
                      
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
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("---------->\t")
                          cat(sprintf(as.character(array_haplotypes_sel)))
                          cat("\n")
                        }
                        
                        Haplotype_PEER_G_HET_haplotypes_sel<-droplevels(Haplotype_PEER_G_HET_haplotypes[which(Haplotype_PEER_G_HET_haplotypes$Haplotype%in%c("HOM_REF",array_haplotypes_sel)),])
                        
                        # if(Condition_DEBUG == 1)
                        # {
                        #   cat("Haplotype_PEER_G_HET_haplotypes_sel\n")
                        #   cat(str(Haplotype_PEER_G_HET_haplotypes_sel))
                        #   cat("\n")
                        # }
                        
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
                        
                        path10<-paste(out,SELECTED_VARS_UPDATED_sel,'/','Haplotypes','/',paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
                        
                        setwd(path10)
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("path10\n")
                          cat(sprintf(as.character(path10)))
                          cat("\n")
                        }
                        
                        if (file.exists(paste("DTU_LM_HET_RESIDUALS_",comparison_2,'.csv',sep='')))
                        {
                          
                          residuals_df = as.data.frame(fread(file=paste("DTU_LM_HET_RESIDUALS_",comparison_2,'.csv',sep=''),sep=",", header=T), stringsAsFactors=F)
                          colnames(residuals_df)[which(colnames(residuals_df) == "ensembl_gene_id")]<-"transcript_id"
                          
                          if(Condition_DEBUG == 1)
                          {
                            cat("residuals_df\n")
                            cat(str(residuals_df))
                            cat("\n")
                          }
                          
                          residuals_df_sel<-residuals_df[which(residuals_df$transcript_id%in%ENST_array),]
                          if(dim(residuals_df_sel)[1] >0)
                          {
                            if(Condition_DEBUG == 1)
                            {
                              cat("residuals_df_sel\n")
                              cat(str(residuals_df_sel))
                              cat("\n")
                            }
                            
                            if(l ==1)
                            {
                              residuals_df_sel.m<-melt(residuals_df_sel, id.vars=c("transcript_id"),value.name = "residuals_Ratio",variable.name = "sample_id")
                              
                              
                            }#l >1
                            else{
                              
                              residuals_df_sel<-residuals_df_sel[,-which(colnames(residuals_df_sel)%in%samples_HOM_REF)]
                              
                              residuals_df_sel.m<-melt(residuals_df_sel, id.vars=c("transcript_id"),value.name = "residuals_Ratio",variable.name = "sample_id")
                              
                              
                            }
                            
                            if(Condition_DEBUG == 1)
                            {
                              cat("residuals_df_sel.m\n")
                              cat(str(residuals_df_sel.m))
                              cat("\n")
                            }
                            
                            
                            
                            list_haplotypes_residuals[[l]]<-residuals_df_sel.m
                            
                          }#dim(residuals_df_sel)[1] >0
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
                        
                        residuals_df_ALL_Haplotypes$transcript_id<-factor(residuals_df_ALL_Haplotypes$transcript_id,
                                                                          levels=rev(levels_transcript_id),
                                                                          ordered=T)
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("residuals_df_ALL_Haplotypes_1\n")
                          cat(str(residuals_df_ALL_Haplotypes))
                          cat("\n")
                          
                          # #####################################################
                          # quit(status = 1)
                        }
                        
                        
                        
                        A<-round(summary(residuals_df_ALL_Haplotypes$residuals_Ratio[!is.na(residuals_df_ALL_Haplotypes$residuals_Ratio)]),2)
                        # A2<-round(summary(residuals_df_ALL_Haplotypes$TPM[!is.na(residuals_df_ALL_Haplotypes$TPM)]),2)
                        
                        residuals_df_ALL_Haplotypes$residuals_Ratio_no_negative<-residuals_df_ALL_Haplotypes$residuals_Ratio+abs(A[1])
                        
                        A3<-round(summary(residuals_df_ALL_Haplotypes$residuals_Ratio_no_negative[!is.na(residuals_df_ALL_Haplotypes$residuals_Ratio_no_negative)]),0)
                        
                        
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("residuals_df_ALL_Haplotypes_2\n")
                          cat(str(residuals_df_ALL_Haplotypes))
                          cat("\n")
                          
                          cat("Summary_residuals_Ratio:\t")
                          cat(sprintf(as.character(names(A))))
                          cat("\n")
                          cat(sprintf(as.character(A)))
                          cat("\n")
                          
                          cat("Summary_residuals_Ratio_no_negative:\t")
                          cat(sprintf(as.character(names(A3))))
                          cat("\n")
                          cat(sprintf(as.character(A3)))
                          cat("\n")
                          
                          
                        }
                        
                        
                        
                        
                        #### Merge with covariates matrix & keep HET
                        
                        
                        
                        residuals_df_ALL_Haplotypes<-merge(residuals_df_ALL_Haplotypes,
                                                           Haplotype_PEER_G_HET_haplotypes,
                                                           by=c("sample_id"))
                        
                        
                        if(Condition_DEBUG == 1)
                        {
                          
                          cat("residuals_df_ALL_Haplotypes_3\n")
                          cat(str(residuals_df_ALL_Haplotypes))
                          cat("\n")
                          
                          # quit(status = 1)
                          
                        }
                        
                        
                        
                        
                        
                        A<-round(summary(residuals_df_ALL_Haplotypes$residuals_Ratio_no_negative[!is.na(residuals_df_ALL_Haplotypes$residuals_Ratio_no_negative)]),2)
                        
                        
                        
                        step<-abs(A[6]-A[1])/10
                        
                        if(step == 0)
                        {
                          
                          step<-0.1
                        }
                        
                        breaks.Rank<-unique(seq(from= A[1], to=A[6]+step,by=step))
                        
                        A<-summary(residuals_df_ALL_Haplotypes$residuals_Ratio_no_negative)
                        
                        
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
                        
                        
                        
                        if(Condition_DEBUG == 1)
                        {
                          cat("labels.Rank2:\t")
                          cat(sprintf(as.character(labels.Rank)))
                          cat("\n")
                          
                          # quit(status = 1)
                        }
                        
                        
                        
                        
                        graph_adjusted_Ratio<-ggplot(data=residuals_df_ALL_Haplotypes,
                                                     aes(x=Haplotype, y=residuals_Ratio_no_negative, fill=Haplotype)) +
                          geom_violin()+
                          stat_summary(fun = median, fun.min = median, fun.max = median,
                                       geom = "crossbar", width = 0.5)+
                          theme_bw()+
                          scale_x_discrete(name=NULL, drop=F)+
                          scale_y_continuous(name="Residuals Ratio Transcript TPM/geneTPM adjusted for covariates",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
                          scale_fill_brewer(palette = "Dark2", drop=F)+
                          ggtitle(paste(unique(residuals_df_ALL_Haplotypes$ensembl_gene_id),unique(residuals_df_ALL_Haplotypes$HGNC),sep=' '))+
                          ggeasy::easy_center_title()
                        
                        # cat("facet_wrap_1:\t")
                        # cat("\n")
                        
                        
                        
                        graph_adjusted_Ratio <- graph_adjusted_Ratio +
                          facet_grid(cols = vars(residuals_df_ALL_Haplotypes$transcript_id), scales = "free") +
                          theme_cowplot(font_size = 14)+
                          theme( strip.background = element_blank(),
                                 strip.placement = "inside",
                                 strip.text = element_text(size=14),
                                 panel.spacing = unit(0.2, "lines"), 
                                 panel.background=element_rect(fill="white"),
                                 panel.border=element_rect(colour="black",size=1),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())+
                          theme(axis.title.y=element_text(size=24, family="sans"),
                                axis.title.x=element_text(size=24, family="sans"),
                                axis.text.y=element_text(angle=0,size=18, color="black", family="sans"),
                                axis.text.x=element_text(angle=45,size=18,vjust=1,hjust=1, color="black", family="sans"),
                                legend.title=element_text(size=16,color="black", family="sans"),
                                legend.text=element_text(size=12,color="black", family="sans"))+
                          theme(legend.position="hidden",legend.title=element_blank(), legend.text = element_text(size=12))
                        
                        #### Print per violin plot gene #2----
                        
                        
                        path7<-paste(out,'FINAL_RESULTS','/','GRAPHICAL_SUMMARY','/', sep='')
                        path9<-paste(path7,SELECTED_VARS_UPDATED_sel,'/', sep='')
                        path9<-paste(path9,'Haplotypes','/', sep='')
                        path10<-paste(path9,paste(SELECTED_VARS_UPDATED_sel,Proxy_array_sel, sep="__"),'/', sep='')
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
                        
                        
                        
                        
                        graph_DEF<-plot_grid(graph_adjusted_Ratio,
                                             nrow = 1,
                                             ncol = 1)
                        
                        
                        
                        
                        svgname<-paste("Violin_plot_DTU_",HGNC_sel,".svg",sep='')
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
                          cat("\n")
                          cat("\n")
                          
                        }
                        
                      }# length(list_haplotypes_residuals) >0
                    }#length(ENST_array) >0
                    
                  }#dim(Results_LogLM_sel_ENSG_sel)[1] >0
                }#z in 1:length(ENSG_array
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
  
 
  graphical_function(opt)
 
 
  
  
  
}


###########################################################################

system.time( main() )
