#!/bin/bash

work_dir=/u/project/zaitlenlab/bballiu/FastGxE/
JobType=1
          # 1: Simulation study
          # 2: Merge GTEx expression files
          # 3: Residualized GTEx expression files for covariates
          # 4: Decompose GTEx expression files and compute PCA
          # 5: Filter duplicated SNPs from genotype file
          # 6: Compute MAF per tissue in GTEx
          # 7: Run MatrixEQTL per tissue in GTEx
          # 8: Get number of SNPs tested per gene and genes tested per SNP in each tissue in GTEx
          # 9: Filter MatrixEQTL results in GTEx for p-value <= threshold
          #10: Run TreeQTL in GTEx
          #11: Summarize eAssociation files across tissues in GTEx
          #12: Prepare files for mashR in GTEx
          #13: Run mashR analysis in GTEx
          #14: eVariants enrichment in genomic features
          #15: Run Alexis Battle method (sn_spMF): make input files
          #16: Run Alexis Battle method (sn_spMF): model selection: grid search for hyperparameters (step 1: narrow down range)
          #17: Run Alexis Battle method (sn_spMF): model selection: grid search for hyperparameters (step 2: refine selection)

# 0: Simulation study
if ([ $JobType -eq 1 ]); then
  nT=5
  I=1e04

  for i in $(seq 1 30); do
    qsub -N sim.$i -o ${work_dir}/logfiles/sim.study.scenario.${i}.o -e ${work_dir}/logfiles/sim.study.scenario.${i}.e -l h_data=16G,h_rt=12:00:00 $work_dir/scripts/submit_job.sh $JobType $nT $I $i $work_dir
  done

fi

# Andrew can you add here the scripts for going from ped to txt expression files? Whenever you get a chance. Just in case we ever have to re-run the analysis.
# if ([ $JobType -eq 2.0 ]); then
# echo "Under construction"
# fi

# Merge expression files
if ([ $JobType -eq 2 ]); then
  qsub -N merge -o ${work_dir}/logfiles/merge.o -e ${work_dir}/logfiles/merge.e -l h_data=64G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType
fi

# Residualized expression files for covariates
if ([ $JobType -eq 3 ]); then
  for i in $(seq 1 49); do
    qsub -N j.$i -o ${work_dir}/logfiles/resexp.${i}.o -e ${work_dir}/logfiles/resexp.${i}.e -l h_data=16G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType $i
  done
fi

# Decompose expression files and compute PCA
if ([ $JobType -eq 4 ]); then
  qsub -N decomp -o ${work_dir}/logfiles/decomp.o -e ${work_dir}/logfiles/decomp.e -l h_data=80G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType
fi

# Andrew can you add here the scripts for processing the vcf files? Maybe merge with the next script that filters the duplicated SNPs. Whenever you get a chance. Just in case we ever have to re-run the analysis.
# if ([ $JobType -eq 5.0 ]); then
# echo "Under construction"
# fi

# Filter duplicated SNPs from genotype file
if ([ $JobType -eq 5 ]); then
    qsub -N remove_dup_snps -o ${work_dir}/logfiles/remove_dup_snps.o -e ${work_dir}/logfiles/remove_dup_snps.e -l h_data=96G,h_rt=4:00:00,highp $work_dir/scripts/submit_job.sh $JobType
fi

# Compute MAF per tissue
if ([ $JobType -eq 6 ]); then
  qsub -N MAF -o ${work_dir}/logfiles/MAF.o -e ${work_dir}/logfiles/MAF.e -l h_data=80G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType
fi

# Run MatrixEQTL by context
if ([ $JobType -eq 7 ]); then
  for exp_scale in $(seq 1 4); do
                    # 1: normalized_expression
                    # 2: normalized_and_residualized_expression
                    # 3: normalized_expression_heterogeneous
                    # 4: normalized_and_residualized_expression_heterogeneous
    for i in $(seq 1 49); do
      qsub -N eQTL.${exp_scale}.${i} -o ${work_dir}/logfiles/eQTL.${exp_scale}.${i}.o -e ${work_dir}/logfiles/eQTL.${exp_scale}.${i}.e -l h_data=64G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType $exp_scale $i
    done
  done

  for exp_scale in $(seq 5 6); do
                       # 5: normalized_expression_homogeneous
                       # 6: normalized_and_residualized_expression_homogeneous
  i=1
      qsub -N eQTL.${exp_scale}.${i} -o ${work_dir}/logfiles/eQTL.${exp_scale}.${i}.o -e ${work_dir}/logfiles/eQTL.${exp_scale}.${i}.e -l h_data=64G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType $exp_scale $i
  done

fi

# Get number of SNPs tested per gene and number of genes tested per SNP in each tissue
if ([ $JobType -eq 8 ]); then

  # Number of SNPs tested per gene for full and het
  nSNPs=1
  cd ${work_dir}/results/eQTL_mapping/MatrixEQTL

  ls -l *.v8.EUR.normalized_and_residualized_expression.all_pairs.txt | awk '{print $9}' > tmp_MatrixEQTL_files

  cat tmp_MatrixEQTL_files | while read MatrixEQTL_file ; do
    tissue_name=$(echo $MatrixEQTL_file | tr "." " " | awk '{print $1}')
    qsub -N nrTests.${tissue_name}.SNPs -o ${work_dir}/logfiles/nrTests.${tissue_name}.SNPs.o -e ${work_dir}/logfiles/nrTests.${tissue_name}.SNPs.e -l h_data=16G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType $tissue_name $nSNPs
  done

  rm tmp_MatrixEQTL_files

  # Number of genes tested per SNPs for full and het
  nSNPs=2
  tissue_name="Whole_Blood"
  qsub -N nrTests.Genes -o ${work_dir}/logfiles/nrTests.Genes.o -e ${work_dir}/logfiles/nrTests.Genes.e -l h_data=16G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType $tissue_name $nSNPs

  # Number of SNPs tested per gene and genes tested per SNPs for hom
  tissue_name="AverageTissue"
  for nSNPs in 3 4; do
    qsub -N nrTests.${tissue_name}.${nSNPs} -o ${work_dir}/logfiles/nrTests.${tissue_name}.${nSNPs}.o -e ${work_dir}/logfiles/nrTests.${tissue_name}.${nSNPs}.e -l h_data=16G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType $tissue_name $nSNPs
  done

fi

# Filter MatrixEQTL results for p-value <= threshold
if ([ $JobType -eq 9 ]); then
  cd ${work_dir}/results/eQTL_mapping/MatrixEQTL
  ls -l *.all_pairs.txt | awk '{print $9}' > tmp_MatrixEQTL_files

  cat tmp_MatrixEQTL_files | while read file ; do
  qsub -N FilterQTLs.${file} -o ${work_dir}/logfiles/FilterQTLs.${file}.o -e ${work_dir}/logfiles/FilterQTLs.${file}.e -l h_data=16G,h_rt=12:00:00,highp $work_dir/scripts/submit_job.sh $JobType ${file}
  done

  rm tmp_MatrixEQTL_files

fi

# Run TreeQTL
if ([ $JobType -eq 10 ]); then
  for exp_scale in $(seq 1 6); do
    qsub -N treeQTL.${exp_scale} -o ${work_dir}/logfiles/treeQTL.${exp_scale}.o -e ${work_dir}/logfiles/treeQTL.${exp_scale}.e -l h_data=80G,h_rt=20:00:00,highp $work_dir/scripts/submit_job.sh $JobType $exp_scale
  done
fi

# Summarize eAssociation files across tissues in GTEx
if ([ $JobType -eq 11 ]); then
  qsub -N sum_tree -o ${work_dir}/logfiles/sum_tree.o -e ${work_dir}/logfiles/sum_tree.e -l h_data=64G,h_rt=3:00:00,highp $work_dir/scripts/submit_job.sh $JobType

fi

# Prepare input files for mashR (Andrew)
if ([ $JobType -eq 12 ]); then
  qsub -N mashR_input -o ${work_dir}/logfiles/mashR_input.o -e ${work_dir}/logfiles/mashR_input.e -l h_data=40G,h_rt=3:00:00,highp $work_dir/scripts/submit_job.sh $JobType

fi

# Run mashR analysis
if ([ $JobType -eq 13 ]); then
  for exp_scale in 2 4; do
    echo Under construction
  done
fi

# eVariants enrichment in genomic features
if ([ $JobType -eq 14 ]); then
  qsub -N enrich -o ${work_dir}/logfiles/enrich.o -e ${work_dir}/logfiles/enrich.e -l h_data=32G,h_rt=15:00:00,highp $work_dir/scripts/submit_job.sh $JobType
fi

# Run Alexis Battle method (sn_spMF): make input files
if ([ $JobType -eq 15 ]); then
  qsub -N sn_spMF_input -o ${work_dir}/logfiles/sn_spMF_input.o -e ${work_dir}/logfiles/sn_spMF_input.e -l h_data=48G,h_rt=20:00:00,highp $work_dir/scripts/submit_job.sh $JobType

fi

# Run Alexis Battle method (sn_spMF): model selection: grid search for hyperparameters (step 1: narrow down range)
if ([ $JobType -eq 16 ]); then

  iterations=20
  for K in 15 20 25
  do
          for alpha1 in 1 10 100 500
          do
                  for lambda1 in 1 10 100 500
                  do
                          hyperparameter_search_run_name=hp.${K}.${alpha1}.${lambda1}
                          qsub -N ${hyperparameter_search_run_name} \
                          -o ${work_dir}/logfiles/${hyperparameter_search_run_name}.o \
                          -e ${work_dir}/logfiles/${hyperparameter_search_run_name}.e \
                          -l h_data=8G,h_rt=95:00:00,highp \
                          $work_dir/scripts/submit_job.sh $JobType ${K} ${alpha1} ${lambda1} ${iterations}

                  done
          done
  done
fi
