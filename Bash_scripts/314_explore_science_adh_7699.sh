#!/bin/bash>

MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4



output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/314_CREATION.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output


#### Rscript

Rscript=/software/R-4.1.0/bin/Rscript


echo "###################################################### explore_tables  #####################################################################"  >> $output
echo "#########################################################################################################################################################################"  >> $output

Rscript_explore_tables=/nfs/users/nfs_m/mt19/Scripts/R/351_Exploration_of_science_adh_7699.R


Table_S3A=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/science_adh_7699/Table_S3A.txt')
Table_S3F=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/science_adh_7699/Table_S3F_REAL.txt')
Table_S3G=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/science_adh_7699/Table_S3G.txt')
Table_S3I=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/science_adh_7699/Table_S3I.txt')
Table_S3K=$(echo '/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/science_adh_7699/Table_S3K.txt')
Table_Supp=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/VARIANT_INTERPRETATION/OVERHAUL_INTERVAL/FINAL_RESULTS/Fig4_pannels/Supp_Table_4_CURATED_Plus_phenotypes.rds")


type=$(echo "explore_tables")


outfile_explore_tables=$(echo "$output_dir""outfile""_""$type"".out")
touch $outfile_explore_tables
echo -n "" > $outfile_explore_tables
name_explore_tables=$(echo "$type""_job")

 
step_mem=$(expr $mem \* 1)
step_pc=$(expr $pc \* 1)

echo "$mem""->""$step_mem"
echo "$pc""->""$step_pc"

 

echo "bsub -G team151 -o $outfile_explore_tables -M $step_mem  -J $name_explore_tables -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$step_mem] rusage[mem=$step_mem] span[hosts=1]\" -n$step_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_explore_tables \\" >> $output
echo "--Table_S3A $Table_S3A \\" >> $output
echo "--Table_S3F $Table_S3F \\" >> $output
echo "--Table_S3G $Table_S3G \\" >> $output
echo "--Table_S3I $Table_S3I \\" >> $output
echo "--Table_S3K $Table_S3K \\" >> $output
echo "--Table_Supp $Table_Supp \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output





bash $output
