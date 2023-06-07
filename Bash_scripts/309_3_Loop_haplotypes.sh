#!/bin/bash>
 
  
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/309_3.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output



MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4

 


# LOOP HAPLOTYPES ###################################################################################################################################################################################
# LOOP HAPLOTYPES ####################################

# LOOP HAPLOTYPES ###################################################################################################################################################################################
# LOOP HAPLOTYPES ###################################################################################################################################################################################
# LOOP HAPLOTYPES ###################################################################################################################################################################################
# LOOP HAPLOTYPES ###################################################################################################################################################################################


SELECTED_VARS=$(echo 'chr10_9032555_C_G,chr1_101744646_C_G,chr1_158613314_G_A,chr1_158638260_G_A,chr4_1008212_C_T,chr14_74218569_C_T,chr12_110331967_G_A,chr2_169313518_C_T,chr12_112092462_G_A,chr12_22610292_G_C,chr12_53309769_A_G,chr9_135857646_C_T,chr7_50444152_G_T,chr14_23587046_A_C,chr14_51132622_T_C,chr6_7143859_C_A,chr14_74630170_A_G,chr1_202129205_G_A,chr8_130641322_C_T,chr15_65174494_A_G,chr15_65174438_C_A,chr16_24761046_G_A,chr16_333719_G_C,chr16_67250992_C_T,chr16_67690688_C_T,chr16_85595360_G_C,chr16_86016328_C_T,chr16_89094897_C_T,chr17_16949211_C_A,chr17_27197056_G_T,chr17_27778073_C_T,chr17_38764524_T_A,chr17_47856909_A_G,chr17_56339594_A_C,chr13_28604007_T_C,chr17_58602131_G_A,chr17_7106378_G_A,chr18_60880701_T_C,chr18_60920854_C_T,chr18_67856078_G_A,chr19_10676941_G_A,chr8_41589736_T_G,chr1_91606142_A_G,chr9_114663385_T_C,chr19_11210157_C_T,chr19_3203962_C_T,chr1_93482787_G_A,chr19_35776481_C_T,chr20_25409287_A_G,chr20_37544151_C_T,chr20_55990370_A_T,chr2_144084356_C_T,chr2_144162105_A_G,chr19_15653669_T_C,chr22_18252442_G_A,chr17_56603493_C_T,chr22_28761148_C_T,chr22_39362450_C_T,chr2_24091099_C_T,chr1_92981236_T_G,chr2_31476771_G_C,chr2_46293826_C_T,chr2_74920648_G_A,chr3_128317978_C_T,chr3_128322617_G_A,chr3_17098399_A_G,chr3_184091102_T_G,chr3_46354444_C_T,chr3_71355240_G_C,chr6_41924998_C_T,chr5_1041433_C_T,chr5_1093511_G_A,chr6_41952511_T_G,chr5_75563535_G_A,chr5_35476470_G_T,chr22_50949811_T_C,chr7_100309180_A_G,chr1_92925654_G_C,chr6_82476412_C_T,chr6_82501719_C_T,chr7_100083971_G_A,chr1_29217311_G_A,chr7_100314474_C_T,chr7_101499930_G_A,chr7_139875230_C_T,chr2_219020958_C_T,chr16_155132_C_T,chr8_130429059_G_A,chr7_99760955_G_C,chr1_198680015_G_A,chr8_90995426_C_T,chr12_111844956_C_T,chr15_64349614_G_A,chr9_135874752_G_A,chr9_135920196_C_T,chr9_20576825_C_A,chr9_4533390_C_G,chr9_5079248_T_C')

queue=$(echo "long")

a=($(echo "$SELECTED_VARS" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        VAR_sel=${i}
#        echo "$VAR_sel"

        VAR_ROUTE=$(echo "$MASTER_ROUTE""/""$VAR_sel""/")


      ### DE LM Haplotypes ###################

      Rscript_DE_LM_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/300_OVERHAUL_DE_partLM_Haplotypes.R


      type=$(echo "DE_LM_Haplotypes""_""$VAR_sel")
      outfile_DE_LM_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_DE_LM_Haplotypes
      echo -n "" > $outfile_DE_LM_Haplotypes
      name_DE_LM_Haplotypes=$(echo "$type""_job")

      GENES_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_GENES_table.txt")
      INTERVAL_GENE_EXP=$(echo "/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/processed_v97/manuel/covariates_peer/ge_matrix_residuals_scaled_final.csv")


      DE_mem=$(expr $mem \* 2)
      DE_pc=$(expr $pc \* 2)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"



      echo "bsub -G team151 -o $outfile_DE_LM_Haplotypes -M $DE_mem -J $name_DE_LM_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_DE_LM_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_table $GENES_table \\" >> $output
      echo "--INTERVAL_GENE_EXP $INTERVAL_GENE_EXP \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output
      
###############       MT_correction_DE_Haplotypes

      Rscript_MT_correction_DE_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/303_OVERHAUL_DE_partMultiple_Testing_Main_Haplotype.R

      type=$(echo "MT_correction_DE_Haplotypes""_""$VAR_sel")
      outfile_MT_correction_DE_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_MT_correction_DE_Haplotypes
      echo -n "" > $outfile_MT_correction_DE_Haplotypes
      name_MT_correction_DE_Haplotypes=$(echo "$type""_job")


      ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
      GENES_PER_BLOCKS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/Allelic_Series_candidates_Allelic_Series_Generation.tsv")
      PCHiC_info=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
      VEP_CSQ=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/VEP_CSQ.csv")

      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


      echo "bsub -G team151 -o $outfile_MT_correction_DE_Haplotypes -M $DE_mem -w\"done($name_DE_LM_Haplotypes)\" -J $name_MT_correction_DE_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#      echo "bsub -G team151 -o $outfile_MT_correction_DE_Haplotypes -M $DE_mem -J $name_MT_correction_DE_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_MT_correction_DE_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
      echo "--PCHiC_info $PCHiC_info \\" >> $output
      echo "--VEP_CSQ $VEP_CSQ \\" >> $output
      echo "--ALL_dB $ALL_dB \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

 
	

 ##############     DTU LM Haplotypes ###################

      Rscript_DTU_logRatio_LM_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/301_OVERHAUL_DTU_part_LogRatio_LM_Haplotypes.R


      type=$(echo "DTU_logRatio_LM_Haplotypes""_""$VAR_sel")
      outfile_DTU_logRatio_LM_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_DTU_logRatio_LM_Haplotypes
      echo -n "" > $outfile_DTU_logRatio_LM_Haplotypes
      name_DTU_logRatio_LM_Haplotypes=$(echo "$type""_job")


      TRANSCRIPTS_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_Transcripts_table.txt")
      INTERVAL_isoform_EXP=$(echo "/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/salmon/tpm_matrix/post_gene_qc/txi_transcript_log2_tmm_tpm_smpls_rmvd_swapped_post_geneqc_01tpm_20perc.csv")
      Transposed_Isoform_Expression=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/Transposed_Isoform_Expression_df.rds")



      DTU_mem=$(expr $mem \* 2)
      DTU_pc=$(expr $pc \* 2)


      echo "$mem""->""$DTU_mem"
      echo "$pc""->""$DTU_pc"


 
      echo "bsub -G hematopoiesis -o $outfile_DTU_logRatio_LM_Haplotypes -M $DTU_mem -J $name_DTU_logRatio_LM_Haplotypes -R\"select[model==Intel_Platinum]\" -R\" select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_DTU_logRatio_LM_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--TRANSCRIPTS_table $TRANSCRIPTS_table \\" >> $output
      echo "--INTERVAL_isoform_EXP $INTERVAL_isoform_EXP \\" >> $output
      echo "--Transposed_Isoform_Expression $Transposed_Isoform_Expression \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

     #   ### MT_correction_DTU_Haplotypes

      Rscript_MT_correction_DTU_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/305_OVERHAUL_DTU_partMultiple_Testing_Haplotype.R

      type=$(echo "MT_correction_DTU_Haplotypes""_""$VAR_sel")
      outfile_MT_correction_DTU_Haplotypes=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_MT_correction_DTU_Haplotypes
      echo -n "" > $outfile_MT_correction_DTU_Haplotypes
        name_MT_correction_DTU_Haplotypes=$(echo "$type""_job")


      ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
      GENES_PER_BLOCKS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/Allelic_Series_candidates_Allelic_Series_Generation.tsv")
      PCHiC_info=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
      VEP_CSQ=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/VEP_CSQ.csv")

      DTU_mem=$(expr $mem \* 1)
      DTU_pc=$(expr $pc \* 1)


      echo "$mem""->""$DTU_mem"
      echo "$pc""->""$DTU_pc"



      echo "bsub -G team151 -o $outfile_MT_correction_DTU_Haplotypes -M $DTU_mem -w\"done($name_DTU_logRatio_LM_Haplotypes)\" -J $name_MT_correction_DTU_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
#      echo "bsub -G team151 -o $outfile_MT_correction_DTU_Haplotypes -M $DTU_mem -w\"done($name_DTU_logRatio_LM_Haplotypes) && done($name_MT_correction_DE_Haplotypes)\" -J $name_MT_correction_DTU_Haplotypes -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_MT_correction_DTU_Haplotypes \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
      echo "--PCHiC_info $PCHiC_info \\" >> $output
      echo "--VEP_CSQ $VEP_CSQ \\" >> $output
      echo "--ALL_dB $ALL_dB \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

      
	

# 	# "chr1_101744646_C_G"
# 	# "chr8_41589736_T_G"
# 	# "chr3_46354444_C_T"
#       # "chr9_135920196_C_T"
#       # "chr3_71355240_G_C"
# #       if [ $VAR_sel == "chr14_51132622_T_C" ]; then

# # #      	  bash $output
# #       	  exit
# #       else

# #       	  echo "Hello_world"" ""$VAR_sel"
#       #       fi


      if [ $VAR_sel == "chr10_9032555_C_G" ]; then
          MT_correction_DTU_Haplotypes_string=$(echo "done($name_MT_correction_DTU_Haplotypes)")
      else
          MT_correction_DTU_Haplotypes_string=$(echo "&& done($name_MT_correction_DTU_Haplotypes)")

      fi


      echo "->>>$MT_correction_DTU_Haplotypes_string"
      arr[${#arr[@]}]="$MT_correction_DTU_Haplotypes_string"


done


done_string=$(echo "\""""${arr[@]}"""\"")
echo "$done_string"
queue=$4

# # #### DE Rscript_PUT_TOGETHER_RESULTS_Haplotypes ----

Rscript_PUT_TOGETHER_RESULTS_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/311_OVERHAUL_DE_PUT_TOGETHER_RESULTS_Haplotypes.R

type=$(echo "PUT_TOGETHER_RESULTS_Haplotypes_DE")
outfile_PUT_TOGETHER_RESULTS_Haplotypes=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_PUT_TOGETHER_RESULTS_Haplotypes
echo -n "" > $outfile_PUT_TOGETHER_RESULTS_Haplotypes
name_PUT_TOGETHER_RESULTS_Haplotypes=$(echo "$type""_job")

DE_mem=$(expr $mem \* 4)
DE_pc=$(expr $pc \* 4)


echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"


echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_Haplotypes -M $DE_mem -w$done_string -J $name_PUT_TOGETHER_RESULTS_Haplotypes -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_Haplotypes -M $DE_mem -J $name_PUT_TOGETHER_RESULTS_Haplotypes -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_PUT_TOGETHER_RESULTS_Haplotypes \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

# #### DE Rscript_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL ----

Rscript_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL=/nfs/users/nfs_m/mt19/Scripts/R/312_OVERHAUL_DE_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL_correction_v2.R

type=$(echo "PUT_TOGETHER_RESULTS_Haplotypes_DE_ALL_BY_ALL")
outfile_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL
echo -n "" > $outfile_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL
name_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL=$(echo "$type""_job")

DE_mem=$(expr $mem \* 8)
DE_pc=$(expr $pc \* 8)


echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"


echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL -M $DE_mem -w$done_string -J $name_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL -M $DE_mem -J $name_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

# #### DTU Rscript_PUT_TOGETHER_RESULTS ----

Rscript_PUT_TOGETHER_RESULTS_Haplotypes_DTU=/nfs/users/nfs_m/mt19/Scripts/R/313_OVERHAUL_DTU_PUT_TOGETHER_RESULTS_Haplotypes_v2.R

type=$(echo "PUT_TOGETHER_RESULTS_Haplotypes_DTU")
outfile_PUT_TOGETHER_RESULTS_Haplotypes_DTU=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_PUT_TOGETHER_RESULTS_Haplotypes_DTU
echo -n "" > $outfile_PUT_TOGETHER_RESULTS_Haplotypes_DTU
name_PUT_TOGETHER_RESULTS_Haplotypes_DTU=$(echo "$type""_job")

DTU_mem=$(expr $mem \* 4)
DTU_pc=$(expr $pc \* 4)


echo "$mem""->""$DTU_mem"
echo "$pc""->""$DTU_pc"


echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_Haplotypes_DTU -M $DTU_mem -w$done_string -J $name_PUT_TOGETHER_RESULTS_Haplotypes_DTU -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_Haplotypes_DTU -M $DTU_mem -J $name_PUT_TOGETHER_RESULTS_Haplotypes_DTU -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_PUT_TOGETHER_RESULTS_Haplotypes_DTU \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

# #### DE graphical summary ----


Rscript_graphical_SUMMARY_DE_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/321_OVERHAUL_DE_Graphical_summary_Haplotypes.R

type=$(echo "graphical_SUMMARY_DE_Haplotypes")
outfile_graphical_SUMMARY_DE_Haplotypes=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_graphical_SUMMARY_DE_Haplotypes
echo -n "" > $outfile_graphical_SUMMARY_DE_Haplotypes
name_graphical_SUMMARY_DE_Haplotypes=$(echo "$type""_job")

DE_mem=$(expr $mem \* 4)
DE_pc=$(expr $pc \* 4)


INTERVAL_GENE_EXP=$(echo "/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/phenotypes/rna_seq/processed_v97/manuel/covariates_peer/ge_matrix_residuals_scaled_final.csv")
Results_INTERVAL_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/ALL_BY_ALL_DE_LM_FPKM_results_Haplotypes.tsv")
Results_BP_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/DE_BP/""ALL_BY_ALL_DE_LM_FPKM_results.tsv")
Supp4_CURATED_WITH_PHENOTYPES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/""Supp_Table_4_CURATED_Plus_phenotypes.rds")

echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"



echo "bsub -G team151 -o $outfile_graphical_SUMMARY_DE_Haplotypes -M $DE_mem -w\"done($name_PUT_TOGETHER_RESULTS_Haplotypes_ALL_BY_ALL)\" -J $name_graphical_SUMMARY_DE_Haplotypes -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_graphical_SUMMARY_DE_Haplotypes -M $DE_mem -J $name_graphical_SUMMARY_DE_Haplotypes -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_graphical_SUMMARY_DE_Haplotypes \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--INTERVAL_GENE_EXP $INTERVAL_GENE_EXP \\" >> $output
echo "--Results_INTERVAL_LM $Results_INTERVAL_LM \\" >> $output
echo "--Results_BP_LM $Results_BP_LM \\" >> $output
echo "--Supp4_CURATED_WITH_PHENOTYPES $Supp4_CURATED_WITH_PHENOTYPES \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



 

############################DTU graphical summary ----



Rscript_graphical_SUMMARY_DTU_Haplotypes=/nfs/users/nfs_m/mt19/Scripts/R/322_OVERHAUL_DTU_Graphical_summary_Haplotypes.R

type=$(echo "graphical_SUMMARY_DTU_Haplotypes")
outfile_graphical_SUMMARY_DTU_Haplotypes=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_graphical_SUMMARY_DTU_Haplotypes
echo -n "" > $outfile_graphical_SUMMARY_DTU_Haplotypes
name_graphical_SUMMARY_DTU_Haplotypes=$(echo "$type""_job")
 
DTU_mem=$(expr $mem \* 8)
DTU_pc=$(expr $pc \* 8)


Results_LogLM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/DTU_LM_logRatio_results_Haplotypes.tsv")
Transposed_Isoform_Expression=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/Transposed_Isoform_Expression_df.rds")
Tappas_gff=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/tappAS/Homo_sapiens_GRCh38_Ensembl_86.gff3")


echo "$mem""->""$DTU_mem"
echo "$pc""->""$DTU_pc"


echo "bsub -G team151 -o $outfile_graphical_SUMMARY_DTU_Haplotypes -M $DTU_mem -w\"done($name_PUT_TOGETHER_RESULTS_Haplotypes_DTU)\" -J $name_graphical_SUMMARY_DTU_Haplotypes -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_graphical_SUMMARY_DTU_Haplotypes -M $DTU_mem -J $name_graphical_SUMMARY_DTU_Haplotypes -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DTU_mem] rusage[mem=$DTU_mem] span[hosts=1]\" -n$DTU_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_graphical_SUMMARY_DTU_Haplotypes \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--Results_LogLM $Results_LogLM \\" >> $output
echo "--Tappas_gff $Tappas_gff \\" >> $output
echo "--Transposed_Isoform_Expression $Transposed_Isoform_Expression \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


bash $output
