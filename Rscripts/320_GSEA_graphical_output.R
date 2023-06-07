

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
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)


graphical_readout = function(option_list)
{
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
  suppressMessages(library("clusterProfiler", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
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
  
  #### READ GSEA_GLOBAL_RESULTS ----
  
  GSEA_GLOBAL_RESULTS = as.data.frame(fread(opt$GSEA_GLOBAL_RESULTS, sep="\t", header=T), stringsAsFactors=F)
  
  cat("GSEA_GLOBAL_RESULTS_0\n")
  cat(str(GSEA_GLOBAL_RESULTS))
  cat("\n")
  
  GSEA_GLOBAL_RESULTS$minus_logpval_adjusted<-round(-1*log10(GSEA_GLOBAL_RESULTS$p.adjust),2)
  
  
  ######################### MASTER LOOP #############################
  
  
  
  SUMMARY_pval_GSEA<-summary(GSEA_GLOBAL_RESULTS$minus_logpval_adjusted)
  
  cat("SUMMARY_minus_logpval_adjusted_GSEA\n")
  cat(sprintf(as.character(names(SUMMARY_pval_GSEA))))
  cat("\n")
  cat(sprintf(as.character(SUMMARY_pval_GSEA)))
  cat("\n")
  
  
  # SELECTED_VARS_UPDATED<-"chr7_101499930_G_A"
  
  
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
    
    ## size scale ----
    
    minuslospval_vector_INTERVAL<-unique(c(Results_INTERVAL_LM$Genome_wide_minuslogpvalue))
    
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
    
    
    breaks.size<-sort(unique(c(1,seq(0,20, by=4))))
    labels.size<-as.character(breaks.size)
    
    cat("labels.size\n")
    cat(str(labels.size))
    cat("\n")
    
    # SELECTED_VARS_UPDATED<-"chr7_101499930_G_A"
    
    
    list_GW<-list()
    
    for(i in 1:length(SELECTED_VARS_UPDATED))
    {
      
      SELECTED_VARS_sel<-SELECTED_VARS_UPDATED[i]
      
      cat("------------------------------------------------------------------------------------------>\t")
      cat(sprintf(as.character(SELECTED_VARS_sel)))
      cat("\n")
      
      if(SELECTED_VARS_sel == "chr2_74920648_G_A")
      {
        Condition_DEBUG <- 1
      }else{

        Condition_DEBUG <- 0
      }
      #### Read files for Residuals violin plots ----
      
      path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
      
      
      setwd(path6)
      
      if(file.exists("DE_RESULTS_FOR_GSEA.rds"))
      {
        GSEA_FILE<-readRDS(file="DE_RESULTS_FOR_GSEA.rds")
        
        if(Condition_DEBUG == 1)
        {
          cat("GSEA_FILE_0\n")
          cat(str(GSEA_FILE))
          cat("\n")
          cat(str(unique(GSEA_FILE$VAR)))
          cat("\n")
          cat(str(unique(GSEA_FILE$ensembl_gene_id)))
          cat("\n")
        }
        
        
        GSEA_GLOBAL_RESULTS_sel<-GSEA_GLOBAL_RESULTS[which(GSEA_GLOBAL_RESULTS$VAR==SELECTED_VARS_sel),]
        
        if(dim(GSEA_GLOBAL_RESULTS_sel)[1] >0)
        {
          GSEA_GLOBAL_RESULTS_sel[order(GSEA_GLOBAL_RESULTS_sel$minus_logpval_adjusted, decreasing=T),]
          
          if(Condition_DEBUG == 1)
          {
            cat("GSEA_GLOBAL_RESULTS_sel_0\n")
            cat(str(GSEA_GLOBAL_RESULTS_sel))
            cat("\n")
            cat(str(unique(GSEA_GLOBAL_RESULTS_sel$VAR)))
            cat("\n")
            cat(str(unique(GSEA_GLOBAL_RESULTS_sel$ID)))
            cat("\n")
          }
          
          GO_ID_array<-unique(GSEA_GLOBAL_RESULTS_sel$ID)
          
          if(Condition_DEBUG == 1)
          {
            cat("GO_ID_array_\n")
            cat(str(GO_ID_array))
            cat("\n")
            
          }
          
          # GO_ID_array<-"GO:0061484"
          
          for(k in 1:length(GO_ID_array))
          {
            GO_ID_array_sel<-GO_ID_array[k]
            
            GSEA_GLOBAL_RESULTS_sel_GO_sel<-GSEA_GLOBAL_RESULTS_sel[which(GSEA_GLOBAL_RESULTS_sel$ID == GO_ID_array_sel),]
            
            minus_logpval_adjusted_GO_sel<-unique(GSEA_GLOBAL_RESULTS_sel_GO_sel$minus_logpval_adjusted)
            
            
            
            if(Condition_DEBUG == 1)
            {
              cat("GSEA_GLOBAL_RESULTS_sel_GO_sel_0\n")
              cat(str(GSEA_GLOBAL_RESULTS_sel_GO_sel))
              cat("\n")
              cat(str(unique(GSEA_GLOBAL_RESULTS_sel_GO_sel$VAR)))
              cat("\n")
              cat(str(unique(GSEA_GLOBAL_RESULTS_sel_GO_sel$ID)))
              cat("\n")
            }
            
            
            Description_sel<-gsub("\\s+","_",unique(GSEA_GLOBAL_RESULTS_sel_GO_sel$Description))
            Description_sel<-gsub(",","_",Description_sel)
            Description_sel<-gsub("-","_",Description_sel)
            Description_sel<-gsub("\\/","_",Description_sel)
            
            GO_ID_array_sel<-gsub(":","_",GO_ID_array_sel)
            
            if(Condition_DEBUG == 1)
            {
              cat("---->\t")
              cat(sprintf(as.character(GO_ID_array_sel)))
              cat("\t")
              cat(sprintf(as.character(Description_sel)))
              cat("\n")
            }
           
            NES_sel<-round(unique(GSEA_GLOBAL_RESULTS_sel_GO_sel$NES),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("NES\t")
              cat(sprintf(as.character(NES_sel)))
              cat("\n")
            }
            
            
            
            ENSG_array<-unlist(strsplit(GSEA_GLOBAL_RESULTS_sel_GO_sel$string_ENSG_core,split=";"))
            
            if(Condition_DEBUG == 1)
            {
              cat("ENSG_array_\n")
              cat(str(ENSG_array))
              cat("\n")
              
            }
            
            GSEA_FILE_ENSG_sel<-GSEA_FILE[which(GSEA_FILE$ensembl_gene_id%in%ENSG_array),]
            GSEA_FILE_ENSG_sel<-GSEA_FILE_ENSG_sel[order(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue, decreasing = T),]
            
            if(Condition_DEBUG == 1)
            {
              cat("GSEA_FILE_ENSG_sel_0\n")
              cat(str(GSEA_FILE_ENSG_sel))
              cat("\n")
              cat(str(unique(GSEA_FILE_ENSG_sel$VAR)))
              cat("\n")
              cat(str(unique(GSEA_FILE_ENSG_sel$ensembl_gene_id)))
              cat("\n")
            }
            
            # quit(status = 1)
            
            ### Significance ----
            
            GSEA_FILE_ENSG_sel$Significance<-"NA"
            
            GSEA_FILE_ENSG_sel$Significance[which(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue >= 1.3)]<-"YES"
            GSEA_FILE_ENSG_sel$Significance[which(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue < 1.3)]<-"NO"
            
            GSEA_FILE_ENSG_sel$Significance<-factor(GSEA_FILE_ENSG_sel$Significance,
                                                    levels=c("NO","YES"),
                                                    ordered=T)
            if(Condition_DEBUG == 1)
            {
              cat("GSEA_FILE_ENSG_sel$Significance\n")
              cat(sprintf(as.character(names(summary(GSEA_FILE_ENSG_sel$Significance)))))
              cat("\n")
              cat(sprintf(as.character(summary(GSEA_FILE_ENSG_sel$Significance))))
              cat("\n")
            }
            
            GSEA_FILE_ENSG_sel_SIG<-GSEA_FILE_ENSG_sel[which(GSEA_FILE_ENSG_sel$Significance == "YES"),]
            
            if(Condition_DEBUG == 1)
            {
              cat("GSEA_FILE_ENSG_sel_SIG_0\n")
              cat(str(GSEA_FILE_ENSG_sel_SIG))
              cat("\n")
              cat(str(unique(GSEA_FILE_ENSG_sel_SIG$VAR)))
              cat("\n")
              cat(str(unique(GSEA_FILE_ENSG_sel_SIG$ensembl_gene_id)))
              cat("\n")
            }
            
            
            if(dim(GSEA_FILE_ENSG_sel_SIG)[1] >4)
            {
              #### RNASeq_source ----
              
              GSEA_FILE_ENSG_sel$RNASeq_source<-"Whole blood"
              
              
              GSEA_FILE_ENSG_sel$RNASeq_source<-factor(GSEA_FILE_ENSG_sel$RNASeq_source,
                                                       levels=c("Whole blood"),
                                                       ordered=T)
              if(Condition_DEBUG == 1)
              {
                cat("GSEA_FILE_ENSG_sel$RNASeq_source\n")
                cat(sprintf(as.character(names(summary(GSEA_FILE_ENSG_sel$RNASeq_source)))))
                cat("\n")
                cat(sprintf(as.character(summary(GSEA_FILE_ENSG_sel$RNASeq_source))))
                cat("\n")
              }
              
              
              ### update break.size
              
              local_minuslogpval_max<-max(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue[!is.na(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue)])
              
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
              
              indx.finite<-is.finite(GSEA_FILE_ENSG_sel$FC)
              
              check.finite<-sum(indx.finite)
              
              if(check.finite >0)
              {
                A<-summary(GSEA_FILE_ENSG_sel$FC[indx.finite])
                
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
                
                breaks.FC_INTERVAL<-sort(unique(round(c(0,candidate_vector,-1,1),1)))
                labels.FC_INTERVAL<-as.character(breaks.FC_INTERVAL)
                
              }else{
                
                breaks.FC_INTERVAL<-sort(unique(round(c(0,-1,1),1)))
                labels.FC_INTERVAL<-as.character(breaks.FC_INTERVAL)
                
              }# length(indx.finite) >0
              
              if(Condition_DEBUG == 1)
              {
                cat("breaks.FC_INTERVAL\n")
                cat(str(breaks.FC_INTERVAL))
                cat("\n")
                
                cat("labels.FC_INTERVAL\n")
                cat(str(labels.FC_INTERVAL))
                cat("\n")
                
                # quit(status = 1)
              }
              
              ### Graph ----
              
              if(dim(GSEA_FILE_ENSG_sel)[1] > 20)
              {
                GSEA_FILE_ENSG_sel<-GSEA_FILE_ENSG_sel[c(1:20),]
              }
              
              if(Condition_DEBUG == 1)
              {
                cat("GSEA_FILE_ENSG_sel_for_dot_plot_REDUCED:\t")
                cat(str(GSEA_FILE_ENSG_sel))
                cat("\n")
                
              }
              
              dotplot_INTERVAL<-ggplot(data=GSEA_FILE_ENSG_sel,
                                       aes(y=HGNC,
                                           x=RNASeq_source)) +
                geom_point(aes(color=Significance,
                               fill=FC,
                               size=Genome_wide_minuslogpvalue), stroke=1, shape=21)+
                scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
                scale_size(range = c(0,20), name='-log10pval',
                           breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
                scale_fill_gradient2(
                  low = "blue", 
                  mid = "white", 
                  high = "red", 
                  midpoint = 0,
                  breaks=breaks.FC_INTERVAL,labels=labels.FC_INTERVAL,
                  limits=c(breaks.FC_INTERVAL[1]-0.01,breaks.FC_INTERVAL[length(breaks.FC_INTERVAL)]+0.01),name=paste('FoldChange','median(HET)/median(HOM_REF)',sep="\n"),na.value = "gray")+
                scale_y_discrete(name=NULL, drop=F)+
                theme_classic()+
                theme(axis.title.y=element_blank(),
                      axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())+
                theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
                      legend.key.height = unit(1.5, 'cm'), #change legend key height
                      legend.key.width = unit(1, 'cm'), #change legend key width
                      legend.title = element_text(size=14), #change legend title font size
                      legend.text = element_text(size=14))+ #change legend text font size
                ggtitle(paste(Description_sel,paste("NES=",NES_sel,sep=''),
                              paste("mlogpval=",minus_logpval_adjusted_GO_sel,sep=''),
                              sep=" "))+
                theme(axis.title=element_text(angle=0,size=18, face="bold", color="black", family="sans"))+
                scale_x_discrete(name=NULL, drop=T)+
                ggeasy::easy_center_title()
              
              # dotplot_INTERVAL<-ggplot(data=GSEA_FILE_ENSG_sel,
              #                          aes(y=HGNC,
              #                              x=RNASeq_source)) +
              #   geom_point(aes(color=Significance,
              #                  fill=FC,
              #                  size=Genome_wide_minuslogpvalue), stroke=1, shape=21)+
              #   scale_color_manual(values=c("gray","black"),name='Significant', drop=F)+
              #   scale_size(range = c(0,20), name='-log10pval',
              #              breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
              #   scale_fill_gradient2(
              #     low = "blue", 
              #     mid = "white", 
              #     high = "red", 
              #     midpoint = 0,
              #     breaks=breaks.FC_INTERVAL,labels=labels.FC_INTERVAL,
              #     limits=c(breaks.FC_INTERVAL[1]-0.01,breaks.FC_INTERVAL[length(breaks.FC_INTERVAL)]+0.01),name=paste('FoldChange','median(HET)/median(HOM_REF)',sep="\n"),na.value = "gray")+
              #   scale_y_discrete(name=NULL, drop=F)+
              #   theme_classic()+
              #   scale_x_discrete(name=NULL, drop=T)
              #   ggeasy::easy_center_title()
              # 
              # dotplot_INTERVAL<-dotplot_INTERVAL+
              #   facet_grid(cols = vars(GSEA_FILE_ENSG_sel$RNASeq_source), scales='free_x', space='free_x', drop=T) +
              #   theme_cowplot(font_size = 14)+
              #   theme( strip.background = element_blank(),
              #          strip.placement = "inside",
              #          strip.text = element_text(size=14),
              #          panel.spacing = unit(0.2, "lines"), 
              #          panel.background=element_rect(fill="white"),
              #          panel.border=element_rect(colour="black",size=1),
              #          panel.grid.major = element_blank(),
              #          panel.grid.minor = element_blank())+
              #   theme(axis.title.y=element_blank(),
              #         axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
              #         axis.text.x=element_blank(),
              #         axis.ticks.x=element_blank())+
              #   theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
              #         legend.key.height = unit(1.5, 'cm'), #change legend key height
              #         legend.key.width = unit(1, 'cm'), #change legend key width
              #         legend.title = element_text(size=14), #change legend title font size
              #         legend.text = element_text(size=14))+ #change legend text font size
              #   ggtitle(paste(Description_sel,paste("NES=",NES_sel,sep=''),
              #                 paste("mlogpval=",minus_logpval_adjusted_GO_sel,sep=''),
              #                 sep=" "))+
              #   theme(axis.title=element_text(angle=0,size=18, face="bold", color="black", family="sans"))
                
              
              graph_DEF<-plot_grid(dotplot_INTERVAL,
                                   nrow = 1,
                                   ncol = 1)
              
              
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
              
              path9<-paste(path8,'GSEA','/', sep='')
              
              if (file.exists(path9)){
                
                
                
                
              } else {
                dir.create(file.path(path9))
                
              }
              
              setwd(path9)
              
              svgname<-paste("GSEA_",minus_logpval_adjusted_GO_sel,"_",Description_sel,"_",GO_ID_array_sel,"_",SELECTED_VARS_sel,".svg",sep='')
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
              
              # #####################################################
              # quit(status = 1)
              
            }#dim(GSEA_FILE_ENSG_sel_SIG)[1] >0
            
            SELECTED_GO_TERMS<-c("GO:0015629","GO:0061484")
            
            FLAG_df<-GSEA_GLOBAL_RESULTS_sel_GO_sel[which(GSEA_GLOBAL_RESULTS_sel_GO_sel$ID%in%SELECTED_GO_TERMS),]
            
            FLAG_GO<-dim(FLAG_df)[1]
            
            if(Condition_DEBUG == 1)
            {
              cat("FLAG_df\n")
              cat(str(FLAG_df))
              cat("\n")
              cat("FLAG_GO\n")
              cat(str(FLAG_GO))
              cat("\n")
             
            }
            
            if(length(FLAG_GO) >0)
            {
              #### RNASeq_source ----
              
              GSEA_FILE_ENSG_sel$RNASeq_source<-"Whole blood"
              
              
              GSEA_FILE_ENSG_sel$RNASeq_source<-factor(GSEA_FILE_ENSG_sel$RNASeq_source,
                                                       levels=c("Whole blood"),
                                                       ordered=T)
              if(Condition_DEBUG == 1)
              {
                cat("GSEA_FILE_ENSG_sel$RNASeq_source\n")
                cat(sprintf(as.character(names(summary(GSEA_FILE_ENSG_sel$RNASeq_source)))))
                cat("\n")
                cat(sprintf(as.character(summary(GSEA_FILE_ENSG_sel$RNASeq_source))))
                cat("\n")
              }
              
              
              ### update break.size
              
              local_minuslogpval_max<-max(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue[!is.na(GSEA_FILE_ENSG_sel$Genome_wide_minuslogpvalue)])
              
              if(Condition_DEBUG == 1)
              {
                cat("local_minuslogpval_max\n")
                cat(str(local_minuslogpval_max))
                cat("\n")
              }
              
              
             
              
              breaks.size_updated<-sort(unique(c(1,seq(0,local_minuslogpval_max, by=1))))
              labels.size_updated<-as.character(breaks.size_updated)
              
              
              
              
              
              ### Use Beta from the full model of logRatio ---
              
              ### Labels beta
              
              indx.finite<-is.finite(GSEA_FILE_ENSG_sel$FC)
              
              check.finite<-sum(indx.finite)
              
              if(check.finite >0)
              {
                A<-summary(GSEA_FILE_ENSG_sel$FC[indx.finite])
                
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
                
                breaks.FC_INTERVAL<-sort(unique(round(c(0,candidate_vector,-1,1),1)))
                labels.FC_INTERVAL<-as.character(breaks.FC_INTERVAL)
                
              }else{
                
                breaks.FC_INTERVAL<-sort(unique(round(c(0,-1,1),1)))
                labels.FC_INTERVAL<-as.character(breaks.FC_INTERVAL)
                
              }# length(indx.finite) >0
              
              if(Condition_DEBUG == 1)
              {
                cat("breaks.FC_INTERVAL\n")
                cat(str(breaks.FC_INTERVAL))
                cat("\n")
                
                cat("labels.FC_INTERVAL\n")
                cat(str(labels.FC_INTERVAL))
                cat("\n")
                
                # quit(status = 1)
              }
              
              ### Graph ----
              
              if(dim(GSEA_FILE_ENSG_sel)[1] > 20)
              {
                GSEA_FILE_ENSG_sel<-GSEA_FILE_ENSG_sel[c(1:20),]
              }
              
              if(Condition_DEBUG == 1)
              {
                cat("GSEA_FILE_ENSG_sel_for_dot_plot_REDUCED:\t")
                cat(str(GSEA_FILE_ENSG_sel))
                cat("\n")
                
              }
              
              dotplot_INTERVAL<-ggplot(data=GSEA_FILE_ENSG_sel,
                                       aes(y=HGNC,
                                           x=RNASeq_source)) +
                geom_point(aes(color=Significance,
                               fill=FC,
                               size=Genome_wide_minuslogpvalue), stroke=1, shape=21)+
                scale_color_manual(values=c("black","black"),name='Significant', drop=F)+
                scale_size(range = c(0,20), name='-log10pval',
                           breaks=breaks.size_updated, labels=labels.size_updated, limits=c(breaks.size_updated[1],breaks.size_updated[length(breaks.size_updated)]))+
                scale_fill_gradient2(
                  low = "blue", 
                  mid = "white", 
                  high = "red", 
                  midpoint = 0,
                  breaks=breaks.FC_INTERVAL,labels=labels.FC_INTERVAL,
                  limits=c(breaks.FC_INTERVAL[1]-0.01,breaks.FC_INTERVAL[length(breaks.FC_INTERVAL)]+0.01),name=paste('FoldChange','median(HET)/median(HOM_REF)',sep="\n"),na.value = "gray")+
                scale_y_discrete(name=NULL, drop=F)+
                theme_classic()+
                theme(axis.title.y=element_blank(),
                      axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())+
                theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
                      legend.key.height = unit(1.5, 'cm'), #change legend key height
                      legend.key.width = unit(1, 'cm'), #change legend key width
                      legend.title = element_text(size=14), #change legend title font size
                      legend.text = element_text(size=14))+ #change legend text font size
                ggtitle(paste(Description_sel,paste("NES=",NES_sel,sep=''),
                              paste("mlogpval=",minus_logpval_adjusted_GO_sel,sep=''),
                              sep=" "))+
                theme(axis.title=element_text(angle=0,size=18, face="bold", color="black", family="sans"))+
                scale_x_discrete(name=NULL, drop=T)+
                ggeasy::easy_center_title()
              
              # dotplot_INTERVAL<-dotplot_INTERVAL+
              #   facet_grid(cols = vars(GSEA_FILE_ENSG_sel$RNASeq_source), scales='free_x', space='free_x', drop=T) +
              #   theme_cowplot(font_size = 14)+
              #   theme( strip.background = element_blank(),
              #          strip.placement = "inside",
              #          strip.text = element_text(size=14),
              #          panel.spacing = unit(0.2, "lines"), 
              #          panel.background=element_rect(fill="white"),
              #          panel.border=element_rect(colour="black",size=1),
              #          panel.grid.major = element_blank(),
              #          panel.grid.minor = element_blank())+
              #   theme(axis.title.y=element_blank(),
              #         axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
              #         axis.text.x=element_blank(),
              #         axis.ticks.x=element_blank())+
              #   theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
              #         legend.key.height = unit(1.5, 'cm'), #change legend key height
              #         legend.key.width = unit(1, 'cm'), #change legend key width
              #         legend.title = element_text(size=14), #change legend title font size
              #         legend.text = element_text(size=14))+ #change legend text font size
              #   ggtitle(paste(Description_sel,paste("NES=",NES_sel,sep=''),
              #                 paste("mlogpval=",minus_logpval_adjusted_GO_sel,sep=''),
              #                 sep=" "))+
              #   theme(axis.title=element_text(angle=0,size=18, face="bold", color="black", family="sans"))
                
              
              graph_DEF<-plot_grid(dotplot_INTERVAL,
                                   nrow = 1,
                                   ncol = 1)
              
              
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
              
              path9<-paste(path8,'GSEA','/', sep='')
              
              if (file.exists(path9)){
                
                
                
                
              } else {
                dir.create(file.path(path9))
                
              }
              
              setwd(path9)
              
              svgname<-paste("GSEA_",minus_logpval_adjusted_GO_sel,"_",Description_sel,"_",GO_ID_array_sel,"_",SELECTED_VARS_sel,".svg",sep='')
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
              
             
              
            }#length(FLAG_GO) >0
            
          }#k in 1:length(GO_ID_array)
          
          # #####################################################
          # quit(status = 1)
          
        }#dim(GSEA_GLOBAL_RESULTS_sel)[1] >0
      }#file.exists("DE_RESULTS_FOR_GSEA.rds")
    }#i in 1:length(SELECTED_VARS)
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
    make_option(c("--GSEA_GLOBAL_RESULTS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Results_INTERVAL_LM"), type="character", default=NULL,
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
  
 
  graphical_readout(opt)

  
}


###########################################################################

system.time( main() )
