#!/bin/bash>


#### Rscript

Rscript=/software/R-4.1.0/bin/Rscript



MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/304_CREATION.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output

### FINAL_CLASSIFICATION ----

Rscript_FINAL_CLASSIFICATION=/nfs/users/nfs_m/mt19/Scripts/R/323_curated_MPRA_HITS_Figures_v5.R

type=$(echo "FINAL_CLASSIFICATION")
outfile_FINAL_CLASSIFICATION=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_FINAL_CLASSIFICATION
echo -n "" > $outfile_FINAL_CLASSIFICATION
name_FINAL_CLASSIFICATION=$(echo "$type""_job")

FINAL_CLASSIF_mem=$(expr $mem \* 2)
FINAL_CLASSIF_pc=$(expr $pc \* 2)



CURATED_TABLE=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/ZZZ_NEW_table_after_NEW_DTU_model.csv")
Results_INTERVAL_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/ALL_BY_ALL_DE_LM_FPKM_results_Main_VARS.tsv")
Results_BP_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/DE_BP/""ALL_BY_ALL_DE_LM_FPKM_results.tsv")
Results_LogLM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/DTU_LM_logRatio_results_Main_VARS.tsv")
Proxy_R2_TABLE=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/Proxy_R2_TABLE_ALL.tsv")
Proxy_CSQ_TABLE=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/Proxy_CSQ_TABLE_ALL.tsv")
CUMMULATIVE_CLASSES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/GLOBAL_ANALYSIS_Feb_2022_FC01/cummulative_plots/cummulative_classification.rds")
genIE_input=$(echo '/nfs/users/nfs_m/mt19/analysis_genIE/input.csv')
genIE_active=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/Table_genIE.tsv")
ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
Prob_Threshold=$(echo "0.1")
AF_and_CSQ_file=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")

echo "$mem""->""$FINAL_CLASSIF_mem"
echo "$pc""->""$FINAL_CLASSIF_pc"



echo "bsub -G team151 -o $outfile_FINAL_CLASSIFICATION -M $FINAL_CLASSIF_mem -J $name_FINAL_CLASSIFICATION -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$FINAL_CLASSIF_mem] rusage[mem=$FINAL_CLASSIF_mem] span[hosts=1]\" -n$FINAL_CLASSIF_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_FINAL_CLASSIFICATION \\" >> $output
echo "--CURATED_TABLE $CURATED_TABLE \\" >> $output
echo "--Results_INTERVAL_LM $Results_INTERVAL_LM \\" >> $output
echo "--Results_BP_LM $Results_BP_LM \\" >> $output
echo "--Results_LogLM $Results_LogLM \\" >> $output
echo "--Proxy_R2_TABLE $Proxy_R2_TABLE \\" >> $output
echo "--Proxy_CSQ_TABLE $Proxy_CSQ_TABLE \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--genIE_input $genIE_input \\" >> $output
echo "--genIE_active $genIE_active \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--Prob_Threshold $Prob_Threshold \\" >> $output
echo "--AF_and_CSQ_file $AF_and_CSQ_file \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output
  
# ### POST_MANUAL_CURATION ----

Rscript_POST_MANUAL_CURATION=/nfs/users/nfs_m/mt19/Scripts/R/324_curated_MPRA_HITS_Figures_v7_POST_MANUAL_CURATION.R

type=$(echo "POST_MANUAL_CURATION")
outfile_POST_MANUAL_CURATION=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_POST_MANUAL_CURATION
echo -n "" > $outfile_POST_MANUAL_CURATION
name_POST_MANUAL_CURATION=$(echo "$type""_job")

FINAL_CLASSIF_mem=$(expr $mem \* 2)
FINAL_CLASSIF_pc=$(expr $pc \* 2)



# ZZZ_in_progress=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/ZZZ_in_progress_reinstated_without_observations.tsv")
# ZZZ_in_progress=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/ZZZ_in_progress_reinstated.tsv")
#ZZZ_in_progress=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/ZZZ_in_progress_reinstated_with_observations.tsv")
#ZZZ_in_progress=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/ZZZ_reclassif_TYMP.tsv")
ZZZ_in_progress=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/ZZZ_table_April_2023.tsv")
Results_INTERVAL_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/ALL_BY_ALL_DE_LM_FPKM_results_Main_VARS.tsv")
Results_BP_LM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/DE_BP/""ALL_BY_ALL_DE_LM_FPKM_results.tsv")
Results_LogLM=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/DTU_LM_logRatio_results_Main_VARS.tsv")
Proxy_R2_TABLE=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/Proxy_R2_TABLE_ALL.tsv")
Proxy_CSQ_TABLE=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/Proxy_CSQ_TABLE_ALL.tsv")
CUMMULATIVE_CLASSES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/GLOBAL_ANALYSIS_Feb_2022_FC01/cummulative_plots/cummulative_classification.rds")
genIE_input=$(echo '/nfs/users/nfs_m/mt19/analysis_genIE/input.csv')
genIE_active=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/Table_genIE.tsv")
ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
Prob_Threshold=$(echo "0.1")
AF_and_CSQ_file=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/V2F_paper/MPRA_E_Plus_ASE_ACTIVE_AT_LEAST_1_genIE_Tier_del_0.1/ER_Labelling_Initial_Selection_with_CSQ_labels.rds")

echo "$mem""->""$FINAL_CLASSIF_mem"
echo "$pc""->""$FINAL_CLASSIF_pc"



echo "bsub -G team151 -o $outfile_POST_MANUAL_CURATION -M $FINAL_CLASSIF_mem -J $name_POST_MANUAL_CURATION -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$FINAL_CLASSIF_mem] rusage[mem=$FINAL_CLASSIF_mem] span[hosts=1]\" -n$FINAL_CLASSIF_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_POST_MANUAL_CURATION \\" >> $output
echo "--ZZZ_in_progress $ZZZ_in_progress \\" >> $output
echo "--Results_INTERVAL_LM $Results_INTERVAL_LM \\" >> $output
echo "--Results_BP_LM $Results_BP_LM \\" >> $output
echo "--Results_LogLM $Results_LogLM \\" >> $output
echo "--Proxy_R2_TABLE $Proxy_R2_TABLE \\" >> $output
echo "--Proxy_CSQ_TABLE $Proxy_CSQ_TABLE \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--genIE_input $genIE_input \\" >> $output
echo "--genIE_active $genIE_active \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--Prob_Threshold $Prob_Threshold \\" >> $output
echo "--AF_and_CSQ_file $AF_and_CSQ_file \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output



### STATS_MECH ----

Rscript_STATS_MECH=/nfs/users/nfs_m/mt19/Scripts/R/325_figure_4_table_v5.R

type=$(echo "STATS_MECH")
outfile_STATS_MECH=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_STATS_MECH
echo -n "" > $outfile_STATS_MECH
name_STATS_MECH=$(echo "$type""_job")

STATS_MECH_mem=$(expr $mem \* 1)
STATS_MECH_pc=$(expr $pc \* 1)


Supp4_Table_CURATED_PLUS_PHENOTYPES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/""Supp_Table_4_CURATED_Plus_phenotypes.rds")
RMV_common=$(echo "chr4_1008212_C_T")
RMV_labels=$(echo 'NOT_SCREENED_MPRA,NOT_SCREENED_genIE,No_RNA_Seq_HET_carriers')
CUMMULATIVE_CLASSES=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/NEW_MPRA/GLOBAL_ANALYSIS_Feb_2022_FC01/cummulative_plots/cummulative_classification.rds")
type=$(echo "Fig4_STATS_MECH")



echo "$mem""->""$STATS_MECH_mem"
echo "$pc""->""$STATS_MECH_pc"



echo "bsub -G team151 -o $outfile_STATS_MECH -M $STATS_MECH_mem -w\"done($name_POST_MANUAL_CURATION)\" -J $name_STATS_MECH -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$STATS_MECH_mem] rusage[mem=$STATS_MECH_mem] span[hosts=1]\" -n$STATS_MECH_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_STATS_MECH -M $STATS_MECH_mem -J $name_STATS_MECH -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$STATS_MECH_mem] rusage[mem=$STATS_MECH_mem] span[hosts=1]\" -n$STATS_MECH_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_STATS_MECH \\" >> $output
echo "--Supp4_Table_CURATED_PLUS_PHENOTYPES $Supp4_Table_CURATED_PLUS_PHENOTYPES \\" >> $output
echo "--RMV_common $RMV_common \\" >> $output
echo "--RMV_labels $RMV_labels \\" >> $output
echo "--CUMMULATIVE_CLASSES $CUMMULATIVE_CLASSES \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


bash $output
