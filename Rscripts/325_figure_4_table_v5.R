
suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("cowplot",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("digest",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("farver",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("labeling",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
library("reshape2",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("ggeasy",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("RColorBrewer",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")

library("viridisLite",lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("cowplot", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggnewscale", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("scales", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

opt = NULL

options(warn=1)

initial_barplot_MECH_vs_CURATION = function(option_list)
{
  
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
  
  
  #### READ SUpp TABLE 4----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  
  
  
  #### RMV DGKQ ----
  
  
  RMV_common = opt$RMV_common
  
  cat("RMV_common_\n")
  cat(sprintf(as.character(RMV_common)))
  cat("\n")
  
  #### RMV labels ----
  
  
  RMV_labels = unlist(strsplit(opt$RMV_labels, split=","))
  
  cat("RMV_labels_\n")
  cat(sprintf(as.character(RMV_labels)))
  cat("\n")
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR == RMV_common),]
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class == RMV_labels[3]),]
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class<-droplevels(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_1:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class))))
  cat("\n")
  
  # #############################################################
  # quit(status = 1)
  
  #### TOTAL Manual_curation Mechanistic_Class
  
  index_TOTAL<-c(which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "Manual_curation"),
                 which(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES) == "Mechanistic_Class"))
  
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES.dt<-data.table(Supp4_Table_CURATED_PLUS_PHENOTYPES, key=c(colnames(Supp4_Table_CURATED_PLUS_PHENOTYPES)[index_TOTAL]))
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES.dt\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES.dt))
  cat("\n")
  
  
  Supp4_Table_Freq_Table_detailed<-as.data.frame(Supp4_Table_CURATED_PLUS_PHENOTYPES.dt[, .(Freq=.N), by=key(Supp4_Table_CURATED_PLUS_PHENOTYPES.dt)],stringsAsFactors=F)
  
  
  cat("Supp4_Table_Freq_Table_detailed\n")
  cat(str(Supp4_Table_Freq_Table_detailed))
  cat("\n")
 
  
  A<-summary(Supp4_Table_Freq_Table_detailed$Freq)
  
  cat("summary(Supp4_Table_Freq_Table_detailed$Freq)\n")
  cat(sprintf(as.character(names(A))))
  cat("\n")
  cat(sprintf(as.character(A)))
  cat("\n")
  
  max_A<-as.numeric(A[6])
  min_A<-as.numeric(A[1])
  
  step_A<-round((max_A-min_A)/10,0)
  step_A<-5
  
  cat("step_A\n")
  cat(str(max_A))
  cat("\n")
  cat(str(min_A))
  cat("\n")
  cat(str(step_A))
  cat("\n")
  
  breaks.Rank<-(seq(0,35,by=5))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
   
 
  
 
  #### Graph
  
  
  
  vector_colors_MPRA<-brewer.pal(3, "Blues")
  
  # Supp4_Table_Freq_Table_detailed %>%
  #   mutate(myaxis = paste0(Manual_curation, "\n", "n=", TOTAL)) %>%
  #   mutate(myaxis=fct_reorder(myaxis,as.numeric(Manuel_Category))) %>%
  
  Supp4_Table_graph<-ggplot(data=Supp4_Table_Freq_Table_detailed,
                     aes(x=Manual_curation, y=Freq)) +
    geom_bar(stat="identity",colour='black', fill="steelblue")+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("# variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    ggeasy::easy_center_title()
  # 
  # Supp4_Table_graph<-Supp4_Table_graph+
  #   theme_cowplot(font_size = 12)+
  #   facet_grid(. ~ Mechanistic_Class, drop=F)+
  #   theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
  #   guides(fill=guide_legend(nrow=3,byrow=TRUE))+
  #   theme(axis.text.x=element_text(angle=90,size=8,hjust=1,color="black", family="sans"))
  
  
  Supp4_Table_graph<-Supp4_Table_graph+
    facet_grid(. ~ Mechanistic_Class, drop=F)+
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
          axis.text.x=element_text(angle=90,size=14,vjust=1, hjust=1, color="black", family="sans"))+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))+ #change legend text font size
    scale_x_discrete(name=NULL, drop=F)
  
  cat("Supp4_Table_graph DONE\n")
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  # cat("path7\n")
  # cat(sprintf(as.character(path7)))
  # cat("\n")
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  setwd(path7)
  
  # setwd(out)
  
  graph_DEF<-plot_grid(Supp4_Table_graph,
                       nrow = 1,
                       ncol = 1)
  
  svglite(paste('Fig4_INITIAL_Mech_class_vs_Manual_curation_classes_BARPLOT','.svg',sep=''), width = 8, height = 8)
  print(graph_DEF)
  dev.off()
  
  # ###############################
  # quit(status = 1)
  # 
}

data_wrangling = function(option_list)
{
  
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
  
  
  #### READ SUpp TABLE 4----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  
 
  
  #### RMV DGKQ ----
  
  
  RMV_common = opt$RMV_common
  
  cat("RMV_common_\n")
  cat(sprintf(as.character(RMV_common)))
  cat("\n")
  
  #### RMV labels ----
  
  
  RMV_labels = unlist(strsplit(opt$RMV_labels, split=","))
  
  cat("RMV_labels_\n")
  cat(sprintf(as.character(RMV_labels)))
  cat("\n")
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR == RMV_common),]
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class == RMV_labels[3]),]
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class<-droplevels(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_1:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class))))
  cat("\n")
  
  
  ##### MPRA CLASS plot -----
  
  MPRA_subset<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$MPRA_CLASS == RMV_labels[1]),]
  
  MPRA_subset$MPRA_CLASS<-droplevels(MPRA_subset$MPRA_CLASS)
  
  cat("MPRA_subset_0:\n")
  cat(str(MPRA_subset))
  cat("\n")
  cat(str(unique(MPRA_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_subset$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_subset$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_subset$MPRA_CLASS))))
  cat("\n")
  
  check<-MPRA_subset[is.na(MPRA_subset$MPRA_CLASS),]
  
  cat("check_0:\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$Mechanistic_Class))))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$MPRA_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$MPRA_CLASS))))
  cat("\n")
  
  
  # ########################################
  # quit(status = 1)
  
  #### TOTAL Manual_curation Mechanistic_Class
  
  index_TOTAL<-c(which(colnames(MPRA_subset) == "Manual_curation"),
                 which(colnames(MPRA_subset) == "Mechanistic_Class"))
  
  
  
  MPRA_subset.dt<-data.table(MPRA_subset, key=c(colnames(MPRA_subset)[index_TOTAL]))
  
  cat("MPRA_subset.dt\n")
  cat(str(MPRA_subset.dt))
  cat("\n")
  
  
  MPRA_Freq_Table_TOTAL<-as.data.frame(MPRA_subset.dt[, .(TOTAL=.N), by=key(MPRA_subset.dt)],stringsAsFactors=F)
  
  
  cat("MPRA_Freq_Table_TOTAL\n")
  cat(str(MPRA_Freq_Table_TOTAL))
  cat("\n")
  
  #### detailed
  
  index_detailed<-c(index_TOTAL,
                 which(colnames(MPRA_subset) == "MPRA_CLASS"))
  
  
  
  MPRA_subset.dt<-data.table(MPRA_subset, key=c(colnames(MPRA_subset)[index_detailed]))
  
  cat("MPRA_subset.dt\n")
  cat(str(MPRA_subset.dt))
  cat("\n")
  
  
  MPRA_Freq_Table_detailed<-as.data.frame(MPRA_subset.dt[, .(Freq=.N), by=key(MPRA_subset.dt)],stringsAsFactors=F)
  
  
  cat("MPRA_Freq_Table_detailed_0\n")
  cat(str(MPRA_Freq_Table_detailed))
  cat("\n")
  
  #### Merge
  
  MPRA_Freq_Table_detailed<-merge(MPRA_Freq_Table_detailed,
                             MPRA_Freq_Table_TOTAL,
                             by=colnames(MPRA_Freq_Table_TOTAL)[which(colnames(MPRA_Freq_Table_TOTAL)%in%colnames(MPRA_Freq_Table_detailed))])
  
  
  cat("MPRA_Freq_Table_detailed_1\n")
  cat(str(MPRA_Freq_Table_detailed))
  cat("\n")
  
  MPRA_Freq_Table_detailed$Perc<-100*(MPRA_Freq_Table_detailed$Freq/MPRA_Freq_Table_detailed$TOTAL)
  
  cat("MPRA_Freq_Table_detailed_2\n")
  cat(str(MPRA_Freq_Table_detailed))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(MPRA_Freq_Table_detailed, file="test.tsv", sep="\t", quote=F,row.names = F)
  
  #### Graph
  
  breaks.Rank<-(seq(0,100,by=10))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  vector_colors_MPRA<-brewer.pal(3, "Blues")
  
  # MPRA_Freq_Table_detailed %>%
  #   mutate(myaxis = paste0(Manual_curation, "\n", "n=", TOTAL)) %>%
  #   mutate(myaxis=fct_reorder(myaxis,as.numeric(Manuel_Category))) %>%
  
  MPRA_graph<-ggplot(data=MPRA_Freq_Table_detailed,
                     aes(x=Manual_curation, y=Perc, fill=MPRA_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("Percentage of variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_fill_manual(values=vector_colors_MPRA,drop=F)+
    ggeasy::easy_center_title()
  
  # MPRA_graph<-MPRA_graph+
  #   theme_cowplot(font_size = 12)+
  #   facet_grid(. ~ Mechanistic_Class, drop=F)+
  #   theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
  #   guides(fill=guide_legend(nrow=3,byrow=TRUE))+
  #   theme(axis.text.x=element_text(angle=90,size=8,hjust=1,color="black", family="sans"))
  
 
  MPRA_graph<-MPRA_graph+
    facet_grid(. ~ Mechanistic_Class, drop=F)+
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
          axis.text.x=element_text(angle=90,size=14,vjust=1, hjust=1, color="black", family="sans"))+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))+ #change legend text font size
    scale_x_discrete(name=NULL, drop=F)
    
  cat("MPRA_graph DONE\n")
  
 
  ##### genIE CLASS plot -----
  
  genIE_subset<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$genIE_CLASS == RMV_labels[2]),]
  
  genIE_subset$genIE_CLASS<-droplevels(genIE_subset$genIE_CLASS)
  
  cat("genIE_subset_0:\n")
  cat(str(genIE_subset))
  cat("\n")
  cat(str(unique(genIE_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(genIE_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(genIE_subset$Mechanistic_Class))))
  cat("\n")
  
  #### TOTAL Manual_curation Mechanistic_Class
  
  index_TOTAL<-c(which(colnames(genIE_subset) == "Manual_curation"),
                 which(colnames(genIE_subset) == "Mechanistic_Class"))
  
  
  
  genIE_subset.dt<-data.table(genIE_subset, key=c(colnames(genIE_subset)[index_TOTAL]))
  
  cat("genIE_subset.dt\n")
  cat(str(genIE_subset.dt))
  cat("\n")
  
  
  genIE_Freq_Table_TOTAL<-as.data.frame(genIE_subset.dt[, .(TOTAL=.N), by=key(genIE_subset.dt)],stringsAsFactors=F)
  
  
  cat("genIE_Freq_Table_TOTAL\n")
  cat(str(genIE_Freq_Table_TOTAL))
  cat("\n")
  
  #### detailed
  
  index_detailed<-c(index_TOTAL,
                    which(colnames(genIE_subset) == "genIE_CLASS"))
  
  
  
  genIE_subset.dt<-data.table(genIE_subset, key=c(colnames(genIE_subset)[index_detailed]))
  
  cat("genIE_subset.dt\n")
  cat(str(genIE_subset.dt))
  cat("\n")
  
  
  genIE_Freq_Table_detailed<-as.data.frame(genIE_subset.dt[, .(Freq=.N), by=key(genIE_subset.dt)],stringsAsFactors=F)
  
  
  cat("genIE_Freq_Table_detailed_0\n")
  cat(str(genIE_Freq_Table_detailed))
  cat("\n")
  
  #### Merge
  
  genIE_Freq_Table_detailed<-merge(genIE_Freq_Table_detailed,
                                  genIE_Freq_Table_TOTAL,
                                  by=colnames(genIE_Freq_Table_TOTAL)[which(colnames(genIE_Freq_Table_TOTAL)%in%colnames(genIE_Freq_Table_detailed))])
  
  
  cat("genIE_Freq_Table_detailed_1\n")
  cat(str(genIE_Freq_Table_detailed))
  cat("\n")
  
  genIE_Freq_Table_detailed$Perc<-100*(genIE_Freq_Table_detailed$Freq/genIE_Freq_Table_detailed$TOTAL)
  
  cat("genIE_Freq_Table_detailed_2\n")
  cat(str(genIE_Freq_Table_detailed))
  cat("\n")
  
  # setwd(out)
  # 
  # write.table(genIE_Freq_Table_detailed, file="test.tsv", sep="\t", quote=F,row.names = F)
  
  #### Graph
  
  breaks.Rank<-(seq(0,100,by=10))
  labels.Rank<-as.character(breaks.Rank)
  
  
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  vector_colors_genIE<-brewer.pal(6, "Paired")[c(5:6)]
  
  
  # genIE_Freq_Table_detailed %>%
  #   mutate(myaxis = paste0(Manual_curation, "\n", "n=", TOTAL)) %>%
  #   mutate(myaxis=fct_reorder(myaxis,as.numeric(Manuel_Category))) %>%
  
  genIE_graph<-ggplot(data=genIE_Freq_Table_detailed,
                     aes(x=Manual_curation, y=Perc, fill=genIE_CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_classic()+
    theme(axis.title.y=element_text(size=18, family="sans"),
          axis.text.y=element_text(angle=0,size=14, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=14,hjust=1,vjust=1, color="black", family="sans"))+
    scale_y_continuous(name=paste("Percentage of variants",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_x_discrete(name=NULL, drop=F)+
    scale_fill_manual(values=vector_colors_genIE,drop=F)+
    ggeasy::easy_center_title()
  
  # genIE_graph<-genIE_graph+
  #   theme_cowplot(font_size = 12)+
  #   facet_grid(. ~ Mechanistic_Class, drop=F)+
  #   theme(legend.position="bottom",legend.title=element_blank(), legend.text = element_text(size=10))+
  #   guides(fill=guide_legend(nrow=3,byrow=TRUE))+
  #   theme(axis.text.x=element_text(angle=90,size=8,hjust=1,color="black", family="sans"))
  
  
  genIE_graph<-genIE_graph+
    facet_grid(. ~ Mechanistic_Class, drop=F)+
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
          axis.text.x=element_text(angle=90,size=14,vjust=1, hjust=1, color="black", family="sans"))+
    theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
          legend.key.height = unit(1.5, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=14))+ #change legend text font size
    scale_x_discrete(name=NULL, drop=F)
  
  cat("genIE_graph DONE\n")
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  # cat("path7\n")
  # cat(sprintf(as.character(path7)))
  # cat("\n")
  
  if (file.exists(path7)){
    
    
    
    
  } else {
    dir.create(file.path(path7))
    
  }
  
  setwd(path7)
  
  # setwd(out)
  
  graph_DEF<-plot_grid(MPRA_graph,
                       nrow = 1,
                       ncol = 1)
  
  svglite(paste('Fig4_Barplot_MPRA_SCREEN_Mech_class_vs_and_Manual_curation_classes','.svg',sep=''), width = 8, height = 8)
  print(graph_DEF)
  dev.off()
  
  cat("MPRA_graph DONE\n")
  
  graph_DEF<-plot_grid(genIE_graph,
                       nrow = 1,
                       ncol = 1)
  
  svglite(paste('Fig4_Barplot_genIE_SCREEN_Mech_class_vs_and_Manual_curation_classes','.svg',sep=''), width = 8, height = 8)
  print(graph_DEF)
  dev.off()
  
  cat("genIE_graph DONE\n")
  
  saveRDS(MPRA_Freq_Table_detailed, file="MPRA_Freq_table.rds")
  saveRDS(genIE_Freq_Table_detailed, file="genIE_Freq_table.rds")
 
}

STATS_MPRA = function(option_list)
{
  # library("RVAideMemoire", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  
  #### READ SUpp TABLE 4----
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)
  
  filename="MPRA_Freq_table.rds"
  File_for_STATS_MPRA<-as.data.frame(readRDS(file=filename), stringsAsFactors=F)
  
  File_for_STATS_MPRA$ID<-row.names(File_for_STATS_MPRA)
  
  
  cat("File_for_STATS_MPRA_0:\n")
  cat(str(File_for_STATS_MPRA))
  cat("\n")
  
  filename="genIE_Freq_table.rds"
  File_for_STATS_genIE<-as.data.frame(readRDS(file=filename), stringsAsFactors=F)
  
  File_for_STATS_genIE$ID<-row.names(File_for_STATS_genIE)
  
  
  cat("File_for_STATS_genIE_0:\n")
  cat(str(File_for_STATS_genIE))
  cat("\n")
  
  
  #### MPRA stats ----
  
  array_levels_MECH_MPRA<-levels(File_for_STATS_MPRA$Mechanistic_Class)
  
  cat("array_levels_MECH_MPRA_0:\n")
  cat(str(array_levels_MECH_MPRA))
  cat("\n")
  
  array_levels_CURATION_MPRA<-levels(File_for_STATS_MPRA$Manual_curation)
  
  cat("array_levels_CURATION_MPRA_0:\n")
  cat(str(array_levels_CURATION_MPRA))
  cat("\n")
  
  MPRA_REST_cat<-c("NO_enhancer_activity","AT_LEAST_1_TILE_with_enhancer_activity")
  MPRA_chosen_cat<-c("AT_LEAST_1_TILE_with_E_Plus_ASE_activity")
  
  Condition_DEBUG <- 0
  
  list_MPRA<-list()
  for(i in 1:length(array_levels_MECH_MPRA))
  {
    MECH_sel<-array_levels_MECH_MPRA[i]
    
    cat("----COMPARISON_BAIT--------->\t")
    cat(sprintf(as.character(MECH_sel)))
    cat("\t")
    
    # REST_MECH<-array_levels_MECH_MPRA[-i]
    
    list_2<-list()
    for(k in 1:length(array_levels_CURATION_MPRA))
    {
      CURATION_sel<-array_levels_CURATION_MPRA[k]
      
      cat("------------->\t")
      cat(sprintf(as.character(CURATION_sel)))
      cat("\n")
      
      # REST_CURATION<-array_levels_CURATION_MPRA[-k]
      
      File_for_STATS_MPRA_sel<-File_for_STATS_MPRA[which(File_for_STATS_MPRA$Mechanistic_Class == MECH_sel &
                                                           File_for_STATS_MPRA$Manual_curation ==   CURATION_sel),]
      if(Condition_DEBUG == 1)
      {
        cat("File_for_STATS_MPRA_sel_0:\n")
        cat(str(File_for_STATS_MPRA_sel))
        cat("\n")
      }
      
      if(dim(File_for_STATS_MPRA_sel)[1] >0)
      {
        MPRA_E_Plus_ASE_sel<-0
        
        FLAG_E_Plus_ASE_sel<-length(File_for_STATS_MPRA_sel$MPRA_CLASS[which(File_for_STATS_MPRA_sel$MPRA_CLASS%in%MPRA_chosen_cat)])
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_E_Plus_ASE_sel:\n")
          cat(str(FLAG_E_Plus_ASE_sel))
          cat("\n")
        }
        
        if(FLAG_E_Plus_ASE_sel >0)
        {
          
          MPRA_E_Plus_ASE_sel<-as.integer(File_for_STATS_MPRA_sel$Freq[which(File_for_STATS_MPRA_sel$MPRA_CLASS%in%MPRA_chosen_cat)])
        }
        
        MPRA_OTHER_LEV_sel<-0
        
        FLAG_OTHER_LEV_sel<-length(File_for_STATS_MPRA_sel$MPRA_CLASS[which(File_for_STATS_MPRA_sel$MPRA_CLASS%in%MPRA_REST_cat)])
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_OTHER_LEV_sel:\n")
          cat(str(FLAG_OTHER_LEV_sel))
          cat("\n")
        }
        
        if(FLAG_OTHER_LEV_sel >0)
        {
          
          MPRA_OTHER_LEV_sel<-as.integer(sum(File_for_STATS_MPRA_sel$Freq[which(File_for_STATS_MPRA_sel$MPRA_CLASS%in%MPRA_REST_cat)]))
        }
        
        
        File_for_STATS_MPRA_NOT_sel<-File_for_STATS_MPRA[-which(File_for_STATS_MPRA$ID%in%File_for_STATS_MPRA_sel$ID),]
        
        # cat("File_for_STATS_MPRA_NOT_sel_0:\n")
        # cat(str(File_for_STATS_MPRA_NOT_sel))
        # cat("\n")
        
        REST_array_levels_MECH_MPRA<-levels(File_for_STATS_MPRA_NOT_sel$Mechanistic_Class)
        REST_array_levels_CURATION_MPRA<-levels(File_for_STATS_MPRA_NOT_sel$Manual_curation)
        
        list_3<-list()
        for(z in 1:length(REST_array_levels_MECH_MPRA))
        {
          REST_MECH_sel<-REST_array_levels_MECH_MPRA[z]
          
          list_4<-list()
          for(v in 1:length(REST_array_levels_CURATION_MPRA))
          {
            REST_CURATION_sel<-REST_array_levels_CURATION_MPRA[v]
            
            File_for_STATS_MPRA_NOT_sel_sel<-File_for_STATS_MPRA_NOT_sel[which(File_for_STATS_MPRA_NOT_sel$Mechanistic_Class == REST_MECH_sel &
                                                                                 File_for_STATS_MPRA_NOT_sel$Manual_curation ==   REST_CURATION_sel),]
            
            if(dim(File_for_STATS_MPRA_NOT_sel_sel)[1] >0)
            {
              cat("-----COMPARISON_PREY-------->\t")
              cat(sprintf(as.character(REST_MECH_sel)))
              cat("\t")
              
              cat("------------->\t")
              cat(sprintf(as.character(REST_CURATION_sel)))
              cat("\n")
              
              if(Condition_DEBUG == 1)
              {
                cat("File_for_STATS_MPRA_NOT_sel_sel_0:\n")
                cat(str(File_for_STATS_MPRA_NOT_sel_sel))
                cat("\n")
              }
              
              MPRA_E_Plus_ASE_NOT_sel_sel<-0
              
              FLAG_E_Plus_ASE_NOT_sel_sel<-length(File_for_STATS_MPRA_NOT_sel_sel$MPRA_CLASS[which(File_for_STATS_MPRA_NOT_sel_sel$MPRA_CLASS%in%MPRA_chosen_cat)])
              
              if(Condition_DEBUG == 1)
              {
                cat("FLAG_E_Plus_ASE_NOT_sel_sel:\n")
                cat(str(FLAG_E_Plus_ASE_NOT_sel_sel))
                cat("\n")
              }
              
              if(FLAG_E_Plus_ASE_NOT_sel_sel >0)
              {
                
                MPRA_E_Plus_ASE_NOT_sel_sel<-as.integer(File_for_STATS_MPRA_NOT_sel_sel$Freq[which(File_for_STATS_MPRA_NOT_sel_sel$MPRA_CLASS%in%MPRA_chosen_cat)])
              }
              
              MPRA_OTHER_LEV_NOT_sel_sel<-0
              
              FLAG_OTHER_LEV_sel<-length(File_for_STATS_MPRA_NOT_sel_sel$MPRA_CLASS[which(File_for_STATS_MPRA_NOT_sel_sel$MPRA_CLASS%in%MPRA_REST_cat)])
              
              if(Condition_DEBUG == 1)
              {
                cat("FLAG_OTHER_LEV_sel:\n")
                cat(str(FLAG_OTHER_LEV_sel))
                cat("\n")
              }
              
              if(FLAG_OTHER_LEV_sel >0)
              {
                
                MPRA_OTHER_LEV_NOT_sel_sel<-as.integer(sum(File_for_STATS_MPRA_NOT_sel_sel$Freq[which(File_for_STATS_MPRA_NOT_sel_sel$MPRA_CLASS%in%MPRA_REST_cat)]))
              }
              
              if(Condition_DEBUG == 1)
              {
                cat("MPRA_E_Plus_ASE_sel:\n")
                cat(str(MPRA_E_Plus_ASE_sel))
                cat("\n")

                cat("MPRA_OTHER_LEV_sel:\n")
                cat(str(MPRA_OTHER_LEV_sel))
                cat("\n")

                cat("MPRA_E_Plus_ASE_NOT_sel_sel:\n")
                cat(str(MPRA_E_Plus_ASE_NOT_sel_sel))
                cat("\n")

                cat("MPRA_OTHER_LEV_NOT_sel_sel:\n")
                cat(str(MPRA_OTHER_LEV_NOT_sel_sel))
                cat("\n")
              }
              
              A.df<-as.data.frame(rbind(cbind(MPRA_OTHER_LEV_sel,MPRA_E_Plus_ASE_sel),
                                        cbind(MPRA_OTHER_LEV_NOT_sel_sel,MPRA_E_Plus_ASE_NOT_sel_sel)))
              
              
              colnames(A.df)<-c("MPRA_inactive","MPRA_active")
              
              tab.chisq.test<-chisq.test(A.df,correct = TRUE)
              
              # cat("tab.chisq.test\n")
              # cat(str(tab.chisq.test))
              # cat("\n")
              pval<-as.numeric(tab.chisq.test$p.value)
              log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
              
              cat("log_pval\n")
              cat(str(log_pval))
              cat("\n")
              
              
              comparison_df<-as.data.frame(cbind(MECH_sel,CURATION_sel,REST_MECH_sel,REST_CURATION_sel,MPRA_OTHER_LEV_sel,MPRA_E_Plus_ASE_sel,MPRA_OTHER_LEV_NOT_sel_sel,MPRA_E_Plus_ASE_NOT_sel_sel,log_pval), stringsAsFactors=F)
              colnames(comparison_df)<-c("MECH_c1","CURATION_c1","MECH_c2","CURATION_c2","MPRA_OTHER_LEV_c1","MPRA_E_Plus_ASE_c1","MPRA_OTHER_LEV_c2","MPRA_E_Plus_ASE_c2","minuslogpval")
              
              if(Condition_DEBUG == 1)
              {
                cat("comparison_df\n")
                cat(str(comparison_df))
                cat("\n")
              }
              
              list_4[[v]]<-comparison_df
              
              # #if(REST_MECH_sel == "Alternative_Transcript_Usage" & REST_CURATION_sel == "Regulated_candidate_effector_gene")
              # #if(REST_MECH_sel == "Alternative_Transcript_Usage" & REST_CURATION_sel == "Regulated_candidate_effector_gene")
              # #if(REST_MECH_sel == "No_regulation" & REST_CURATION_sel == "Unknown")
              # #if(REST_MECH_sel == "No_regulation" & REST_CURATION_sel == "Proxy_Non_CODING_with_regulated_candidate_gene")
              # if(REST_MECH_sel == "Transcriptional" & REST_CURATION_sel == "Unknown")
              # {
                # quit(status = 1)
              #   
              # }
              
            }#dim(File_for_STATS_MPRA_NOT_sel_sel)[1] >0
          }#v in 1:length(REST_array_levels_CURATION_MPRA)
          
          if(length(list_4) >0)
          {
            df_v <- as.data.frame(data.table::rbindlist(list_4, fill=T), stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              cat("df_v\n")
              cat(str(df_v))
              cat("\n")
            }
            
            list_3[[z]]<-df_v 
          }
        }#z in 1:length(REST_array_levels_MECH_MPRA)
        
        if(length(list_3) >0)
        {
          df_z <- as.data.frame(data.table::rbindlist(list_3, fill=T), stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("df_z\n")
            cat(str(df_z))
            cat("\n")
          }
          list_2[[k]]<-df_z
        }
      }#dim(File_for_STATS_MPRA_sel)[1] >0
    }#k in 1:length(array_levels_CURATION_MPRA)
    
    if(length(list_2) >0)
    {
      df_k <- as.data.frame(data.table::rbindlist(list_2, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("df_k\n")
        cat(str(df_k))
        cat("\n")
      }
      list_MPRA[[i]]<-df_k
    }
    
  }#i in 1:length(array_levels_MECH_MPRA)
  
  if(length(list_MPRA) >0)
  {
    df_MPRA <- as.data.frame(data.table::rbindlist(list_MPRA, fill=T), stringsAsFactors=F)
    
    cat("df_MPRA\n")
    cat(str(df_MPRA))
    cat("\n")
    
    path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
    
    setwd(path7)
    
    write.table(df_MPRA,file="CHiSq_MPRA_class.tsv",sep="\t",quote=F,row.names =F)
    
  }#length(list_MPRA) >0
}

STATS_genIE = function(option_list)
{
  # library("RVAideMemoire", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
  
  
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
  
  
  #### READ SUpp TABLE 4----
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)
  
  filename="genIE_Freq_table.rds"
  File_for_STATS_genIE<-as.data.frame(readRDS(file=filename), stringsAsFactors=F)
  
  File_for_STATS_genIE$ID<-row.names(File_for_STATS_genIE)
  
  
  cat("File_for_STATS_genIE_0:\n")
  cat(str(File_for_STATS_genIE))
  cat("\n")
  
  filename="genIE_Freq_table.rds"
  File_for_STATS_genIE<-as.data.frame(readRDS(file=filename), stringsAsFactors=F)
  
  File_for_STATS_genIE$ID<-row.names(File_for_STATS_genIE)
  
  
  cat("File_for_STATS_genIE_0:\n")
  cat(str(File_for_STATS_genIE))
  cat("\n")
  
  
  #### genIE stats ----
  
  array_levels_MECH_genIE<-levels(File_for_STATS_genIE$Mechanistic_Class)
  
  cat("array_levels_MECH_genIE_0:\n")
  cat(str(array_levels_MECH_genIE))
  cat("\n")
  
  array_levels_CURATION_genIE<-levels(File_for_STATS_genIE$Manual_curation)
  
  cat("array_levels_CURATION_genIE_0:\n")
  cat(str(array_levels_CURATION_genIE))
  cat("\n")
  
  genIE_REST_cat<-c("genIE_INACTIVE")
  genIE_chosen_cat<-c("genIE_ACTIVE")
  
  Condition_DEBUG <- 0
  
  list_genIE<-list()
  for(i in 1:length(array_levels_MECH_genIE))
  {
    MECH_sel<-array_levels_MECH_genIE[i]
    
    cat("----COMPARISON_BAIT--------->\t")
    cat(sprintf(as.character(MECH_sel)))
    cat("\t")
    
    # REST_MECH<-array_levels_MECH_genIE[-i]
    
    list_2<-list()
    for(k in 1:length(array_levels_CURATION_genIE))
    {
      CURATION_sel<-array_levels_CURATION_genIE[k]
      
      cat("------------->\t")
      cat(sprintf(as.character(CURATION_sel)))
      cat("\n")
      
      # REST_CURATION<-array_levels_CURATION_genIE[-k]
      
      File_for_STATS_genIE_sel<-File_for_STATS_genIE[which(File_for_STATS_genIE$Mechanistic_Class == MECH_sel &
                                                           File_for_STATS_genIE$Manual_curation ==   CURATION_sel),]
      if(Condition_DEBUG == 1)
      {
        cat("File_for_STATS_genIE_sel_0:\n")
        cat(str(File_for_STATS_genIE_sel))
        cat("\n")
      }
      
      if(dim(File_for_STATS_genIE_sel)[1] >0)
      {
        genIE_ACTIVE_sel<-0
        
        FLAG_ACTIVE_sel<-length(File_for_STATS_genIE_sel$genIE_CLASS[which(File_for_STATS_genIE_sel$genIE_CLASS%in%genIE_chosen_cat)])
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_ACTIVE_sel:\n")
          cat(str(FLAG_ACTIVE_sel))
          cat("\n")
        }
        
        if(FLAG_ACTIVE_sel >0)
        {
          
          genIE_ACTIVE_sel<-as.integer(File_for_STATS_genIE_sel$Freq[which(File_for_STATS_genIE_sel$genIE_CLASS%in%genIE_chosen_cat)])
        }
        
        genIE_OTHER_LEV_sel<-0
        
        FLAG_OTHER_LEV_sel<-length(File_for_STATS_genIE_sel$genIE_CLASS[which(File_for_STATS_genIE_sel$genIE_CLASS%in%genIE_REST_cat)])
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_OTHER_LEV_sel:\n")
          cat(str(FLAG_OTHER_LEV_sel))
          cat("\n")
        }
        
        if(FLAG_OTHER_LEV_sel >0)
        {
          
          genIE_OTHER_LEV_sel<-as.integer(sum(File_for_STATS_genIE_sel$Freq[which(File_for_STATS_genIE_sel$genIE_CLASS%in%genIE_REST_cat)]))
        }
        
        
        File_for_STATS_genIE_NOT_sel<-File_for_STATS_genIE[-which(File_for_STATS_genIE$ID%in%File_for_STATS_genIE_sel$ID),]
        
        # cat("File_for_STATS_genIE_NOT_sel_0:\n")
        # cat(str(File_for_STATS_genIE_NOT_sel))
        # cat("\n")
        
        REST_array_levels_MECH_genIE<-levels(File_for_STATS_genIE_NOT_sel$Mechanistic_Class)
        REST_array_levels_CURATION_genIE<-levels(File_for_STATS_genIE_NOT_sel$Manual_curation)
        
        list_3<-list()
        for(z in 1:length(REST_array_levels_MECH_genIE))
        {
          REST_MECH_sel<-REST_array_levels_MECH_genIE[z]
          
          list_4<-list()
          for(v in 1:length(REST_array_levels_CURATION_genIE))
          {
            REST_CURATION_sel<-REST_array_levels_CURATION_genIE[v]
            
            File_for_STATS_genIE_NOT_sel_sel<-File_for_STATS_genIE_NOT_sel[which(File_for_STATS_genIE_NOT_sel$Mechanistic_Class == REST_MECH_sel &
                                                                                 File_for_STATS_genIE_NOT_sel$Manual_curation ==   REST_CURATION_sel),]
            
            if(dim(File_for_STATS_genIE_NOT_sel_sel)[1] >0)
            {
              cat("-----COMPARISON_PREY-------->\t")
              cat(sprintf(as.character(REST_MECH_sel)))
              cat("\t")
              
              cat("------------->\t")
              cat(sprintf(as.character(REST_CURATION_sel)))
              cat("\n")
              
              if(Condition_DEBUG == 1)
              {
                cat("File_for_STATS_genIE_NOT_sel_sel_0:\n")
                cat(str(File_for_STATS_genIE_NOT_sel_sel))
                cat("\n")
              }
              
              genIE_ACTIVE_NOT_sel_sel<-0
              
              FLAG_ACTIVE_NOT_sel_sel<-length(File_for_STATS_genIE_NOT_sel_sel$genIE_CLASS[which(File_for_STATS_genIE_NOT_sel_sel$genIE_CLASS%in%genIE_chosen_cat)])
              
              if(Condition_DEBUG == 1)
              {
                cat("FLAG_ACTIVE_NOT_sel_sel:\n")
                cat(str(FLAG_ACTIVE_NOT_sel_sel))
                cat("\n")
              }
              
              if(FLAG_ACTIVE_NOT_sel_sel >0)
              {
                
                genIE_ACTIVE_NOT_sel_sel<-as.integer(File_for_STATS_genIE_NOT_sel_sel$Freq[which(File_for_STATS_genIE_NOT_sel_sel$genIE_CLASS%in%genIE_chosen_cat)])
              }
              
              genIE_OTHER_LEV_NOT_sel_sel<-0
              
              FLAG_OTHER_LEV_sel<-length(File_for_STATS_genIE_NOT_sel_sel$genIE_CLASS[which(File_for_STATS_genIE_NOT_sel_sel$genIE_CLASS%in%genIE_REST_cat)])
              
              if(Condition_DEBUG == 1)
              {
                cat("FLAG_OTHER_LEV_sel:\n")
                cat(str(FLAG_OTHER_LEV_sel))
                cat("\n")
              }
              
              if(FLAG_OTHER_LEV_sel >0)
              {
                
                genIE_OTHER_LEV_NOT_sel_sel<-as.integer(sum(File_for_STATS_genIE_NOT_sel_sel$Freq[which(File_for_STATS_genIE_NOT_sel_sel$genIE_CLASS%in%genIE_REST_cat)]))
              }
              
              if(Condition_DEBUG == 1)
              {
                cat("genIE_ACTIVE_sel:\n")
                cat(str(genIE_ACTIVE_sel))
                cat("\n")
                
                cat("genIE_OTHER_LEV_sel:\n")
                cat(str(genIE_OTHER_LEV_sel))
                cat("\n")
                
                cat("genIE_ACTIVE_NOT_sel_sel:\n")
                cat(str(genIE_ACTIVE_NOT_sel_sel))
                cat("\n")
                
                cat("genIE_OTHER_LEV_NOT_sel_sel:\n")
                cat(str(genIE_OTHER_LEV_NOT_sel_sel))
                cat("\n")
              }
              
              A.df<-as.data.frame(rbind(cbind(genIE_OTHER_LEV_sel,genIE_ACTIVE_sel),
                                        cbind(genIE_OTHER_LEV_NOT_sel_sel,genIE_ACTIVE_NOT_sel_sel)))
              
              
              colnames(A.df)<-c("genIE_inactive","genIE_active")
              
              tab.chisq.test<-chisq.test(A.df,correct = TRUE)
              
              # cat("tab.chisq.test\n")
              # cat(str(tab.chisq.test))
              # cat("\n")
              pval<-as.numeric(tab.chisq.test$p.value)
              log_pval<-as.numeric(round(-1*log10(tab.chisq.test$p.value),2))
              
              cat("log_pval\n")
              cat(str(log_pval))
              cat("\n")
              
              
              comparison_df<-as.data.frame(cbind(MECH_sel,CURATION_sel,REST_MECH_sel,REST_CURATION_sel,genIE_OTHER_LEV_sel,genIE_ACTIVE_sel,genIE_OTHER_LEV_NOT_sel_sel,genIE_ACTIVE_NOT_sel_sel,log_pval), stringsAsFactors=F)
              colnames(comparison_df)<-c("MECH_c1","CURATION_c1","MECH_c2","CURATION_c2","genIE_OTHER_LEV_c1","genIE_ACTIVE_c1","genIE_OTHER_LEV_c2","genIE_ACTIVE_c2","minuslogpval")
              
              if(Condition_DEBUG == 1)
              {
                cat("comparison_df\n")
                cat(str(comparison_df))
                cat("\n")
              }
              
              list_4[[v]]<-comparison_df
              
              # #if(REST_MECH_sel == "Alternative_Transcript_Usage" & REST_CURATION_sel == "Regulated_candidate_effector_gene")
              # #if(REST_MECH_sel == "Alternative_Transcript_Usage" & REST_CURATION_sel == "Regulated_candidate_effector_gene")
              # #if(REST_MECH_sel == "No_regulation" & REST_CURATION_sel == "Unknown")
              # #if(REST_MECH_sel == "No_regulation" & REST_CURATION_sel == "Proxy_Non_CODING_with_regulated_candidate_gene")
              # if(REST_MECH_sel == "Transcriptional" & REST_CURATION_sel == "Unknown")
              # {
              # quit(status = 1)
              #   
              # }
              
            }#dim(File_for_STATS_genIE_NOT_sel_sel)[1] >0
          }#v in 1:length(REST_array_levels_CURATION_genIE)
          
          if(length(list_4) >0)
          {
            df_v <- as.data.frame(data.table::rbindlist(list_4, fill=T), stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              cat("df_v\n")
              cat(str(df_v))
              cat("\n")
            }
            
            list_3[[z]]<-df_v 
          }
        }#z in 1:length(REST_array_levels_MECH_genIE)
        
        if(length(list_3) >0)
        {
          df_z <- as.data.frame(data.table::rbindlist(list_3, fill=T), stringsAsFactors=F)
          
          if(Condition_DEBUG == 1)
          {
            cat("df_z\n")
            cat(str(df_z))
            cat("\n")
          }
          list_2[[k]]<-df_z
        }
      }#dim(File_for_STATS_genIE_sel)[1] >0
    }#k in 1:length(array_levels_CURATION_genIE)
    
    if(length(list_2) >0)
    {
      df_k <- as.data.frame(data.table::rbindlist(list_2, fill=T), stringsAsFactors=F)
      
      if(Condition_DEBUG == 1)
      {
        cat("df_k\n")
        cat(str(df_k))
        cat("\n")
      }
      list_genIE[[i]]<-df_k
    }
    
  }#i in 1:length(array_levels_MECH_genIE)
  
  if(length(list_genIE) >0)
  {
    df_genIE <- as.data.frame(data.table::rbindlist(list_genIE, fill=T), stringsAsFactors=F)
    
    cat("df_genIE\n")
    cat(str(df_genIE))
    cat("\n")
    
    path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
    
    setwd(path7)
    
    write.table(df_genIE,file="CHiSq_genIE_class.tsv",sep="\t",quote=F,row.names =F)
    
  }#length(list_genIE) >0
}

cummulative_PLOTS_E_Plus_ASE_and_Log_Rank_test = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
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
  
  
  #### READ SUpp TABLE 4----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  
  
  
  #### RMV DGKQ ----
  
  
  RMV_common = opt$RMV_common
  
  cat("RMV_common_\n")
  cat(sprintf(as.character(RMV_common)))
  cat("\n")
  
  #### RMV labels ----
  
  
  RMV_labels = unlist(strsplit(opt$RMV_labels, split=","))
  
  cat("RMV_labels_\n")
  cat(sprintf(as.character(RMV_labels)))
  cat("\n")
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR == RMV_common),]
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class == RMV_labels[3]),]
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class<-droplevels(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_1:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class))))
  cat("\n")
  
  
  ##### MPRA CLASS plot -----
  
  MPRA_subset<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$MPRA_CLASS == RMV_labels[1]),]
  
  MPRA_subset$MPRA_CLASS<-droplevels(MPRA_subset$MPRA_CLASS)
  
  cat("MPRA_subset_0:\n")
  cat(str(MPRA_subset))
  cat("\n")
  cat(str(unique(MPRA_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_subset$Mechanistic_Class))))
  cat("\n")
  
  
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  
  cat("CUMMULATIVE_CLASSES_0:\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
  cat("\n")
  
  CUMMULATIVE_CLASSES$comparison_VAR<-gsub("^chr","",CUMMULATIVE_CLASSES$VAR)
  
  Condition_DEBUG <- 1
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES\n")
    cat(str(CUMMULATIVE_CLASSES))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted<-CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$carried_variants == CUMMULATIVE_CLASSES$comparison_VAR),]
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted<-droplevels(CUMMULATIVE_CLASSES_restricted[which(CUMMULATIVE_CLASSES_restricted$Cell_Type == "ALL_CT"),])
  
  CUMMULATIVE_CLASSES_restricted$enhancer_classif<-"NA"
  CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif<-"NA"
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_1\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
    cat("\n")
  }
  
  
  ########### Build df ---------------
  
  n_MPRA_TILES<-5
  VAR_vector<-as.character(MPRA_subset$VAR)
  Mech_vector<-as.character(MPRA_subset$Mechanistic_Class)
  Manual_vector<-as.character(MPRA_subset$Manual_curation)
  
  n_VAR<-length(VAR_vector)
  
  
  # Log_rank_df<-data.frame(matrix(ncol=6,nrow=n_MPRA_TILES*n_VAR, 
  #                                        dimnames=list(NULL, c("VAR", "Mechanistic_Class",
  #                                                              "Manual_curation","TILE",
  #                                                              "ACTIVE","GROUP"))),
  #                                 stringsAsFactors = F)
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(Mech_vector,n_MPRA_TILES),
                                   rep(Manual_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","Mechanistic_Class","Manual_curation","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
  }
  
  Log_rank_df$GROUP<-paste(Log_rank_df$Mechanistic_Class,Log_rank_df$Manual_curation, sep="|")
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$GROUP))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$GROUP)))))
    cat("\n")
  }
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_PRE\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[which(CUMMULATIVE_CLASSES_restricted$VAR%in%Log_rank_df$VAR),c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                                                                       which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_POST\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by="VAR")
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$E_Plus_ASE_CLASS_TILES))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$E_Plus_ASE_CLASS_TILES)))))
    cat("\n")
  }
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
  }
  
  
  #### Log Rank test ----
  
  Condition_DEBUG <- 0
  
  list_DEF<-list()
  
  for(i in 1 :length(unique(Log_rank_df$GROUP)))
  {
    BAIT<-unique(Log_rank_df$GROUP)[i]
    
    cat("----------->\t")
    cat(sprintf(as.character(BAIT)))
    cat("\t")
    
    list_PREY<-list()
    
    for(k in 1 :length(unique(Log_rank_df$GROUP)))
    {
    
      PREY<-unique(Log_rank_df$GROUP)[k]
      
      cat("----------->\t")
      cat(sprintf(as.character(PREY)))
      cat("\n")
      
      Mini_comparison<-c(BAIT,PREY)
      
      Log_rank_df_subset<-Log_rank_df[which(Log_rank_df$GROUP%in%Mini_comparison),]
      
      FLAG_diversity<-length(unique(Log_rank_df_subset$GROUP))
      
      if(Condition_DEBUG == 1)
      {
        cat("FLAG_diversity_0\n")
        cat(str(FLAG_diversity))
        cat("\n")
      }
      
      if(FLAG_diversity >1)
      {
        Mini_comparison_c1<-as.character(Mini_comparison[1])
        Mini_comparison_c2<-as.character(Mini_comparison[2])
        
        FLAG_diversity_2<-length(unique(Log_rank_df_subset$ACTIVE))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_2_0\n")
          cat(str(FLAG_diversity_2))
          cat("\n")
        }
        
        if(FLAG_diversity_2 >1)
        {
          if(Condition_DEBUG == 1)
          {
            cat("----------->\t")
            cat(sprintf(as.character(Mini_comparison_c1)))
            cat("\t")
            cat(sprintf(as.character(Mini_comparison_c2)))
            cat("\n")
            
            cat("Log_rank_df_subset_0\n")
            cat(str(Log_rank_df_subset))
            cat("\n")
            cat(str(unique(Log_rank_df_subset$VAR)))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_subset$ACTIVE))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(Log_rank_df_subset$ACTIVE)))))
            cat("\n")
          }
          
          #, "less", "greater"
          
          STATS_test<-logrank.test(
            Log_rank_df_subset$TILE,
            Log_rank_df_subset$ACTIVE,
            Log_rank_df_subset$GROUP,
            alternative = c("two.sided"),
            rho = 0,
            gamma = 0,
            event_time_weights = NULL
          )
          
          # if(Condition_DEBUG == 1)
          # {
          #   cat("STATS_test_0\n")
          #   cat(str(STATS_test))
          #   cat("\n")
          # }
          
          
          p_value<-STATS_test$test$p
          minus_log_val<-round(-1*log10(p_value),2)
          
          if(Condition_DEBUG == 1)
          {
            cat("minus_log_val_0\n")
            cat(str(minus_log_val))
            cat("\n")
          }
          
          
          a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
          
          
          colnames(a.df)<-c("GROUP_c1","GROUP_c2","minus_logpval")
          a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
          
          if(Condition_DEBUG == 1)
          {
            cat("a.df\n")
            cat(str(a.df))
            cat("\n")
          }
          
          list_PREY[[k]]<-a.df
          
        }#FLAG_diversity_2 >1
      }#FLAG_diversity >1
    }# k in 1 :length(unique(Log_rank_df$GROUP)
    
    Condition_DEBUG <- 0
    
    if(length(list_PREY) >0)
    {
      
      PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
      
      if(Condition_DEBUG == 1)
      {
        cat("PREY_df_0\n")
        cat(str(PREY_df))
        cat("\n")
        #quit(status = 1)
      }
      
      list_DEF[[i]]<-PREY_df
      
    }# length(list_PREY) >0
    
    Condition_DEBUG <- 0
    
   
  }#i in 1 :length(unique(Log_rank_df$GROUP))
  
  
  Condition_DEBUG <- 1
  
  
  FINAL_df = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
  
  if(Condition_DEBUG == 1)
  {
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    #quit(status = 1)
  }
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)

  write.table(FINAL_df, file="MPRA_log_Rank_test_STAT.tsv", sep="\t", quote=F,row.names = F)
  
  
  #### GRAPH GROUP_sel ----
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_REMEMBER\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$GROUP))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$GROUP)))))
    cat("\n")
    #quit(status = 1)
  }
  
  GROUP_sel<-c('Transcriptional_Regulation|R_in_candidate','No_Regulation_Detected|Unknown_mechanism','Alternative_Transcript_Usage|R_in_candidate','Transcriptional_Regulation|R_in_non_candidate','No_Regulation_Detected|Regulatory_proxy')#,'No_Regulation_Detected|Coding_proxy')
  
  vector_colors_MPRA<-brewer.pal(length(GROUP_sel), "Set1")
  
  Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$GROUP%in%GROUP_sel),]
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_sel_0\n")
    cat(str(Log_rank_df_sel))
    cat("\n")
    # quit(status = 1)
  }
  
  
  Log_rank_df_sel.dt<-data.table(Log_rank_df_sel, key=c("TILE","GROUP"))
  
  Freq_TOTAL<-as.data.frame(Log_rank_df_sel.dt[,.(TOTAL=.N),by=key(Log_rank_df_sel.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-as.data.frame(Log_rank_df_sel.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df_sel.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(Freq_TOTAL,
                    Freq_table,
                    by=c("TILE","GROUP"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$TOTAL),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  
  
  #### ACTUAL GRAPH ---- 
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)
  
  v_parameter<-"NA" 
  
  indx_GROUP<-which(colnames(Freq_table) == "GROUP")
  indx.yaxis<-which(colnames(Freq_table) == "Perc")
  mycols <- vector_colors_MPRA[c(5,2,4,3,1,6)]
  
  pdfname<-paste("Cummulative_frequency_","E_Plus_ASE_activity_","GROUP_selection",".pdf",sep='')
  pdf(file=pdfname, width=5, height=4, pointsize=12)
  
  par(mai=c(0.9,0.9,0.3,0.2))
  lab <- as.character(unique(Freq_table[,indx_GROUP]))
  
  cat("lab\n")
  cat(sprintf(as.character(lab)))
  cat("\n")
  
  
  plot(Freq_table$TILE, Freq_table[,indx.yaxis],
       ty="n", xlab="TILES with E_Plus_ASE activity",
       ylab="Cummulative % of variants MPRA Active in at least 1 Cell Type",
       axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,70))
  
  
  
  # cat("Hello_world1\n")
  
  # points(Freq_table$TILE, Freq_table[,indx.yaxis], col="darkgrey", pch=19)
  
  # cat("Hello_world2\n")
  # 
  # 
  for (i in 1:length(lab))  {
    
    cat("Hello_world3\t")
    cat(sprintf(as.character(lab[i])))
    cat("\n")
    
    ind <- which(Freq_table[,indx_GROUP]==lab[i])
    
    cat("------------------->ind\n")
    cat(str(ind))
    cat("\n")
    color_sel<-mycols[i]
    
    cat(sprintf(as.character(color_sel)))
    cat("\n")
    
    points(Freq_table$TILE[ind], Freq_table[,indx.yaxis][ind],
           pch=19, col=color_sel)
    lines(Freq_table$TILE[ind], Freq_table[,indx.yaxis][ind],
          lty=1,lwd=3, col=color_sel)
    
    cat("END_3\n")
    
    # points(Freq_table$TILE[ind], Freq_table[,indx.yaxis][ind], pch=19, col=color_sel)
  }#i in 1:length(lab)
  
  cat("END_4\n")
  
  
  legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
  axis(1, at=seq(0,70))
  axis(2, las=1)
  
  dev.off()
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)
  
  write.table(Freq_table, file="MPRA_Cummulative_Freq_table.tsv", sep="\t", quote=F,row.names = F)
}












cummulative_PLOTS_E_Plus_ASE_and_Log_Rank_test_Multilineage_vs_lineage_restricted = function(option_list)
{
  suppressMessages(library("nph", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
  
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
  
  
  #### READ SUpp TABLE 4----
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-readRDS(file=opt$Supp4_Table_CURATED_PLUS_PHENOTYPES)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_0:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  
  
  
  
  #### RMV DGKQ ----
  
  
  RMV_common = opt$RMV_common
  
  cat("RMV_common_\n")
  cat(sprintf(as.character(RMV_common)))
  cat("\n")
  
  #### RMV labels ----
  
  
  RMV_labels = unlist(strsplit(opt$RMV_labels, split=","))
  
  cat("RMV_labels_\n")
  cat(sprintf(as.character(RMV_labels)))
  cat("\n")
  
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR == RMV_common),]
  Supp4_Table_CURATED_PLUS_PHENOTYPES<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class == RMV_labels[3]),]
  
  Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class<-droplevels(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)
  
  
  cat("Supp4_Table_CURATED_PLUS_PHENOTYPES_1:\n")
  cat(str(Supp4_Table_CURATED_PLUS_PHENOTYPES))
  cat("\n")
  cat(str(unique(Supp4_Table_CURATED_PLUS_PHENOTYPES$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(Supp4_Table_CURATED_PLUS_PHENOTYPES$Mechanistic_Class))))
  cat("\n")
  
  
  ##### MPRA CLASS plot -----
  
  MPRA_subset<-Supp4_Table_CURATED_PLUS_PHENOTYPES[-which(Supp4_Table_CURATED_PLUS_PHENOTYPES$MPRA_CLASS == RMV_labels[1]),]
  
  MPRA_subset$MPRA_CLASS<-droplevels(MPRA_subset$MPRA_CLASS)
  
  cat("MPRA_subset_0:\n")
  cat(str(MPRA_subset))
  cat("\n")
  cat(str(unique(MPRA_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(MPRA_subset$Mechanistic_Class)))))
  cat("\n")
  cat(sprintf(as.character(summary(MPRA_subset$Mechanistic_Class))))
  cat("\n")
  
  
  ### CUMMULATIVE_CLASSES #----
  
  CUMMULATIVE_CLASSES<-readRDS(file=opt$CUMMULATIVE_CLASSES)
  
  
  cat("CUMMULATIVE_CLASSES_0:\n")
  cat(str(CUMMULATIVE_CLASSES))
  cat("\n")
  cat(str(unique(CUMMULATIVE_CLASSES$VAR)))
  cat("\n")
  
  CUMMULATIVE_CLASSES$comparison_VAR<-gsub("^chr","",CUMMULATIVE_CLASSES$VAR)
  
  Condition_DEBUG <- 1
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES\n")
    cat(str(CUMMULATIVE_CLASSES))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted<-CUMMULATIVE_CLASSES[which(CUMMULATIVE_CLASSES$carried_variants == CUMMULATIVE_CLASSES$comparison_VAR),]
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_0\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
  }
  
  
  
  
  CUMMULATIVE_CLASSES_restricted<-droplevels(CUMMULATIVE_CLASSES_restricted[which(CUMMULATIVE_CLASSES_restricted$Cell_Type != "ALL_CT"),])
  
  CUMMULATIVE_CLASSES_restricted$enhancer_classif<-"NA"
  CUMMULATIVE_CLASSES_restricted$E_Plus_ASE_classif<-"NA"
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_1\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(CUMMULATIVE_CLASSES_restricted$Cell_Type)))))
    cat("\n")
    cat(sprintf(as.character(summary(CUMMULATIVE_CLASSES_restricted$Cell_Type))))
    cat("\n")
  }
  
  
  
  ########### Build df ---------------
  
  Cell_Type_levels<-levels(CUMMULATIVE_CLASSES_restricted$Cell_Type)
  n_Cell_Type_levels<-length(levels(CUMMULATIVE_CLASSES_restricted$Cell_Type))
  n_MPRA_TILES<-5
  VAR_vector<-as.character(MPRA_subset$VAR)
  Mech_vector<-as.character(MPRA_subset$Mechanistic_Class)
  Manual_vector<-as.character(MPRA_subset$Manual_curation)
  Multi_Lineage_vector<-as.character(MPRA_subset$Multi_Lineage)
  
  
  n_VAR<-length(VAR_vector)
  
  if(Condition_DEBUG == 1)
  {
    cat("Parameters\n")
    cat(str(Cell_Type_levels))
    cat("\n")
    cat(str(n_Cell_Type_levels))
    cat("\n")
    cat(str(n_MPRA_TILES))
    cat("\n")
    cat(str(VAR_vector))
    cat("\n")
    cat(str(Mech_vector))
    cat("\n")
    cat(str(Manual_vector))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Multi_Lineage_vector))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Multi_Lineage_vector)))))
    cat("\n")
  }
  
  # ##############################################
  # quit(status = 1)
  
  # Log_rank_df<-data.frame(matrix(ncol=6,nrow=n_MPRA_TILES*n_VAR, 
  #                                        dimnames=list(NULL, c("VAR", "Mechanistic_Class",
  #                                                              "Manual_curation","TILE",
  #                                                              "ACTIVE","GROUP"))),
  #                                 stringsAsFactors = F)
  
  
  
  Log_rank_df<-as.data.frame(cbind(rep(VAR_vector,n_MPRA_TILES),
                                   rep(Mech_vector,n_MPRA_TILES),
                                   rep(Manual_vector,n_MPRA_TILES),
                                   rep(Multi_Lineage_vector,n_MPRA_TILES),
                                   c(rep(1,n_VAR),rep(2,n_VAR),rep(3,n_VAR),rep(4,n_VAR),rep(5,n_VAR)),
                                   rep("NA",n_MPRA_TILES*n_VAR)),
                             stringsAsFactors=F)
  
  colnames(Log_rank_df)<-c("VAR","Mechanistic_Class","Manual_curation","Multi_Lineage","TILE","ACTIVE")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_0\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
   
   
  }
  
  
  
  
 
  
  Log_rank_df$GROUP<-paste(Log_rank_df$Mechanistic_Class,Log_rank_df$Manual_curation,Log_rank_df$Multi_Lineage, sep="|")
  Log_rank_df$TILE<-as.integer(Log_rank_df$TILE)
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_1\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$GROUP))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$GROUP)))))
    cat("\n")
  }
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_PRE\n")
    cat(str(CUMMULATIVE_CLASSES_restricted))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted$VAR)))
    cat("\n")
  }
  
  CUMMULATIVE_CLASSES_restricted_subset<-unique(CUMMULATIVE_CLASSES_restricted[which(CUMMULATIVE_CLASSES_restricted$VAR%in%Log_rank_df$VAR),c(which(colnames(CUMMULATIVE_CLASSES_restricted) == "VAR"),
                                                                                                                                              which(colnames(CUMMULATIVE_CLASSES_restricted) == "Cell_Type"),
                                                                                                                                              which(colnames(CUMMULATIVE_CLASSES_restricted) == "E_Plus_ASE_CLASS_TILES"))])
  
  if(Condition_DEBUG == 1)
  {
    cat("CUMMULATIVE_CLASSES_restricted_subset_POST\n")
    cat(str(CUMMULATIVE_CLASSES_restricted_subset))
    cat("\n")
    cat(str(unique(CUMMULATIVE_CLASSES_restricted_subset$VAR)))
    cat("\n")
  }
  
  Log_rank_df<-merge(Log_rank_df,
                     CUMMULATIVE_CLASSES_restricted_subset,
                     by="VAR")
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$E_Plus_ASE_CLASS_TILES))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$E_Plus_ASE_CLASS_TILES)))))
    cat("\n")
    
    
  }
  
 
  
  ################# Classification of TILES ----------------
  
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE <= Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"1"
  Log_rank_df$ACTIVE[which(Log_rank_df$TILE > Log_rank_df$E_Plus_ASE_CLASS_TILES)]<-"0"
  
  Log_rank_df$ACTIVE<-as.numeric(Log_rank_df$ACTIVE)
  
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_2\n")
    cat(str(Log_rank_df))
    cat("\n")
    cat(str(unique(Log_rank_df$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Log_rank_df$ACTIVE))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Log_rank_df$ACTIVE)))))
    cat("\n")
    
    setwd(out)
    write.table(file="test.tsv", Log_rank_df,sep="\t",quote=F, row.names = F)
  }
  
 
  
  
  #### Log Rank test ----
  
  Condition_DEBUG <- 1
  
  
  CT_array<-levels(Log_rank_df$Cell_Type)
  
  cat("CT_array_0\n")
  cat(str(CT_array))
  cat("\n")
  
  list_DEF_CT<-list()
  
  for(iteration_CT_array in 1:length(CT_array))
  {
    CT_array_sel<-CT_array[iteration_CT_array]
    
    cat("----------->\t")
    cat(sprintf(as.character(CT_array_sel)))
    cat("\t")
    
    Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$Cell_Type == CT_array_sel),]
    
    if(Condition_DEBUG == 1)
    {
      cat("Log_rank_df_sel_0\n")
      cat(str(Log_rank_df_sel))
      cat("\n")
    }
    
    list_DEF<-list()
    
    for(i in 1 :length(unique(Log_rank_df_sel$GROUP)))
    {
      BAIT<-unique(Log_rank_df_sel$GROUP)[i]
      
      # cat("----------->\t")
      cat(sprintf(as.character(BAIT)))
      cat("\t")
      # ##########################
      # quit(status = 1)
      
      list_PREY<-list()
      
      for(k in 1 :length(unique(Log_rank_df_sel$GROUP)))
      {
        
        PREY<-unique(Log_rank_df_sel$GROUP)[k]
        
        cat("----------->\t")
        cat(sprintf(as.character(PREY)))
        cat("\n")
        
        Mini_comparison<-c(BAIT,PREY)
        
        Log_rank_df_sel_subset<-Log_rank_df_sel[which(Log_rank_df_sel$GROUP%in%Mini_comparison),]
        
        FLAG_diversity<-length(unique(Log_rank_df_sel_subset$GROUP))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_diversity_0\n")
          cat(str(FLAG_diversity))
          cat("\n")
        }
        
        if(FLAG_diversity >1)
        {
          Mini_comparison_c1<-as.character(Mini_comparison[1])
          Mini_comparison_c2<-as.character(Mini_comparison[2])
          
          FLAG_diversity_2<-length(unique(Log_rank_df_sel_subset$ACTIVE))
          
          if(Condition_DEBUG == 1)
          {
            cat("FLAG_diversity_2_0\n")
            cat(str(FLAG_diversity_2))
            cat("\n")
          }
          
          if(FLAG_diversity_2 >1)
          {
            if(Condition_DEBUG == 1)
            {
              cat("----------->\t")
              cat(sprintf(as.character(Mini_comparison_c1)))
              cat("\t")
              cat(sprintf(as.character(Mini_comparison_c2)))
              cat("\n")
              
              cat("Log_rank_df_sel_subset_0\n")
              cat(str(Log_rank_df_sel_subset))
              cat("\n")
              cat(str(unique(Log_rank_df_sel_subset$VAR)))
              cat("\n")
              cat(sprintf(as.character(names(summary(as.factor(Log_rank_df_sel_subset$ACTIVE))))))
              cat("\n")
              cat(sprintf(as.character(summary(as.factor(Log_rank_df_sel_subset$ACTIVE)))))
              cat("\n")
            }
            
            #, "less", "greater"
            
            STATS_test<-logrank.test(
              Log_rank_df_sel_subset$TILE,
              Log_rank_df_sel_subset$ACTIVE,
              Log_rank_df_sel_subset$GROUP,
              alternative = c("two.sided"),
              rho = 0,
              gamma = 0,
              event_time_weights = NULL
            )
            
            # if(Condition_DEBUG == 1)
            # {
            #   cat("STATS_test_0\n")
            #   cat(str(STATS_test))
            #   cat("\n")
            # }
            
            
            p_value<-STATS_test$test$p
            minus_log_val<-round(-1*log10(p_value),2)
            
            if(Condition_DEBUG == 1)
            {
              cat("minus_log_val_0\n")
              cat(str(minus_log_val))
              cat("\n")
            }
            
            
            a.df<-as.data.frame(cbind(Mini_comparison_c1,Mini_comparison_c2,minus_log_val), stringsAsFactors=F)
            
            
            colnames(a.df)<-c("GROUP_c1","GROUP_c2","minus_logpval")
            a.df$minus_logpval<-as.numeric(a.df$minus_logpval)
            
            if(Condition_DEBUG == 1)
            {
              cat("a.df\n")
              cat(str(a.df))
              cat("\n")
            }
            
            list_PREY[[k]]<-a.df
            
          }#FLAG_diversity_2 >1
        }#FLAG_diversity >1
      }# k in 1 :length(unique(Log_rank_df_sel$GROUP)
      
      Condition_DEBUG <- 0
      
      if(length(list_PREY) >0)
      {
        
        PREY_df = unique(as.data.frame(data.table::rbindlist(list_PREY, fill=T), stringsAsFactors=F))
        
        if(Condition_DEBUG == 1)
        {
          cat("PREY_df_0\n")
          cat(str(PREY_df))
          cat("\n")
          #quit(status = 1)
        }
        
        list_DEF[[i]]<-PREY_df
        
      }# length(list_PREY) >0
      
      Condition_DEBUG <- 0
      
      
    }#i in 1 :length(unique(Log_rank_df_sel$GROUP))
    
    
    Condition_DEBUG <- 1
    
    
    FINAL_df = unique(as.data.frame(data.table::rbindlist(list_DEF, fill=T), stringsAsFactors=F))
    FINAL_df$Cell_Type<-CT_array_sel
    
    if(Condition_DEBUG == 1)
    {
      cat("FINAL_df_0\n")
      cat(str(FINAL_df))
      cat("\n")
      #quit(status = 1)
    }
    
    Condition_DEBUG <- 0
    
    list_DEF_CT[[iteration_CT_array]]<-FINAL_df
  }#iteration_CT_array
  
  Condition_DEBUG <- 1
  
  FINAL_df_CT = unique(as.data.frame(data.table::rbindlist(list_DEF_CT, fill=T), stringsAsFactors=F))
  
  if(Condition_DEBUG == 1)
  {
    cat("FINAL_df_CT_0\n")
    cat(str(FINAL_df_CT))
    cat("\n")
    #quit(status = 1)
  }
  
 
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)
  
  write.table(FINAL_df, file="MPRA_log_Rank_test_STAT_Lineage_edition_per_CT.tsv", sep="\t", quote=F,row.names = F)
  
  #quit(status = 1)
  
  #### GRAPH GROUP_sel ----
  
  #GROUP_sel<-c('Transcriptional_Regulation|R_in_candidate','No_Regulation_Detected|Unknown_mechanism','Alternative_Transcript_Usage|R_in_candidate','Transcriptional_Regulation|R_in_non_candidate','No_Regulation_Detected|Coding_proxy')
  
  
  GROUP_sel<-c('Transcriptional_Regulation|R_in_candidate|Lineage_restricted',
               'Transcriptional_Regulation|R_in_candidate|Multi_Lineage',
               'No_Regulation_Detected|Unknown_mechanism|Lineage_restricted',
               'No_Regulation_Detected|Unknown_mechanism|Multi_Lineage',
               'Alternative_Transcript_Usage|Coding_proxy|Multi_Lineage')
  
  vector_colors_MPRA<-brewer.pal(length(GROUP_sel), "Set1")
  
  Log_rank_df_sel<-Log_rank_df[which(Log_rank_df$GROUP%in%GROUP_sel),]
  
  if(Condition_DEBUG == 1)
  {
    cat("Log_rank_df_sel_0\n")
    cat(str(Log_rank_df_sel))
    cat("\n")
    #quit(status = 1)
  }
  
  
  Log_rank_df_sel.dt<-data.table(Log_rank_df_sel, key=c("TILE","GROUP","Cell_Type"))
  
  Freq_TOTAL<-as.data.frame(Log_rank_df_sel.dt[,.(TOTAL=.N),by=key(Log_rank_df_sel.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_TOTAL_0\n")
    cat(str(Freq_TOTAL))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-as.data.frame(Log_rank_df_sel.dt[,.(Freq=sum(ACTIVE)),by=key(Log_rank_df_sel.dt)], stringsAsFactors=F)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_0\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table<-merge(Freq_TOTAL,
                    Freq_table,
                    by=c("TILE","GROUP","Cell_Type"),
                    all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_1\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  Freq_table$Perc<-round(100*(Freq_table$Freq/Freq_table$TOTAL),2)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("Freq_table_2\n")
    cat(str(Freq_table))
    cat("\n")
    #quit(status = 1)
  }
  
  path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
  
  
  setwd(path7)
  
  write.table(Freq_table, file=paste("MPRA_Cummulative_Freq_table_","ALL_Cell_Types",".tsv",sep=''), sep="\t", quote=F,row.names = F)
  
  
  #### ACTUAL GRAPH ---- 
  CT_array<-levels(Log_rank_df$Cell_Type)
  
  cat("CT_array_0\n")
  cat(str(CT_array))
  cat("\n")
  
  for(iteration_CT_array in 1:length(CT_array))
  {
    
    CT_array_sel<-CT_array[iteration_CT_array]
    
    cat("----------->\t")
    cat(sprintf(as.character(CT_array_sel)))
    cat("\t")
    
    Freq_table_sel<-Freq_table[which(Freq_table$Cell_Type == CT_array_sel),]
    
    path7<-paste(out,'FINAL_RESULTS','/','Fig4_pannels','/',type,'/', sep='')
    
    
    setwd(path7)
    
    v_parameter<-"NA" 
    
    indx_GROUP<-which(colnames(Freq_table_sel) == "GROUP")
    indx.yaxis<-which(colnames(Freq_table_sel) == "Perc")
    mycols <- vector_colors_MPRA[c(5,2,4,3,1,6)]
    
    pdfname<-paste("Cummulative_frequency_","E_Plus_ASE_activity_","GROUP_selection_",CT_array_sel,".pdf",sep='')
    pdf(file=pdfname, width=5, height=4, pointsize=12)
    
    par(mai=c(0.9,0.9,0.3,0.2))
    lab <- as.character(unique(Freq_table_sel[,indx_GROUP]))
    
    cat("lab\n")
    cat(sprintf(as.character(lab)))
    cat("\n")
    
    
    plot(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis],
         ty="n", xlab="TILES with E_Plus_ASE activity",
         ylab="Cummulative % of variants MPRA Active in at least 1 Cell Type",
         axes=F, cex.lab=1.2, cex.lab=1.3, ylim=c(0,70))
    
    
    
    # cat("Hello_world1\n")
    
    # points(Freq_table_sel$TILE, Freq_table_sel[,indx.yaxis], col="darkgrey", pch=19)
    
    # cat("Hello_world2\n")
    # 
    # 
    for (i in 1:length(lab))  {
      
      cat("Hello_world3\t")
      cat(sprintf(as.character(lab[i])))
      cat("\n")
      
      ind <- which(Freq_table_sel[,indx_GROUP]==lab[i])
      
      cat("------------------->ind\n")
      cat(str(ind))
      cat("\n")
      color_sel<-mycols[i]
      
      cat(sprintf(as.character(color_sel)))
      cat("\n")
      
      points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
             pch=19, col=color_sel)
      lines(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind],
            lty=1,lwd=3, col=color_sel)
      
      cat("END_3\n")
      
      # points(Freq_table_sel$TILE[ind], Freq_table_sel[,indx.yaxis][ind], pch=19, col=color_sel)
    }#i in 1:length(lab)
    
    cat("END_4\n")
    
    
    legend("topright", legend=lab, fill=mycols, border=mycols, bty="n")
    axis(1, at=seq(0,70))
    axis(2, las=1)
    
    dev.off()
    
   
    
  }#iteration_CT_array in 1:length(CT_array)
  
  
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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RMV_common"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--RMV_labels"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUMMULATIVE_CLASSES"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Supp4_Table_CURATED_PLUS_PHENOTYPES"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "137_MPRA_normalization_and_filtering_Rscript_v2.R
                        --regular_table FILE.txt
                        --replicas charac
                        --type1 type1
                        --type2 type2
                        --pvalThreshold integer
                        --FDRThreshold integer
                        --EquivalenceTable FILE.txt
                        --sharpr2Threshold charac",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  
  initial_barplot_MECH_vs_CURATION(opt)
  data_wrangling(opt)
  STATS_MPRA(opt)
  STATS_genIE(opt)
  cummulative_PLOTS_E_Plus_ASE_and_Log_Rank_test(opt)
  cummulative_PLOTS_E_Plus_ASE_and_Log_Rank_test_Multilineage_vs_lineage_restricted(opt)
  
}
  
  
  
 

###########################################################################

system.time( main() )
