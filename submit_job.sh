#!/bin/bash

echo "Loading dependencies"

######## BRUNA ########
echo "Loading dependencies"
. /u/local/Modules/default/init/modules.sh
module load R/3.6.1

######## ANDREW ########
# . /u/home/a/andrewlu/project-zaitlenlab/tools/conda_AL/anaconda3/etc/profile.d/conda.sh
# conda activate r-env-fastgxe

######## ANDREW RUNNING ALEXIS BATTLE METHOD ########
# . /u/local/Modules/default/init/modules.sh
# module load R/3.5.1

JobType=$1

# Simulation study
if ([ $JobType -eq 1 ]); then

nT=$2
I=$3
i=$4
work_dir=$5
echo Scenario $i with $nT tissues and $I iterations

my_script=$work_dir/scripts/simulation_study.R

R --vanilla --slave -f $my_script --args $nT $I $i $work_dir
fi

# Merge expression files
if ([ $JobType -eq 2 ]); then
my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/merge_expression_files.R
R --vanilla --slave -f $my_script
fi

# Residualized expression files
if ([ $JobType -eq 3 ]); then
  i=$2
  echo Running tissue $i
  my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/residualize_expression_for_covariates.R
  R --vanilla --slave -f $my_script --args $i
fi

# Decompose expression files and compute PCA
if ([ $JobType -eq 4 ]); then
  my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/decompose_expression.R
  R --vanilla --slave -f $my_script
fi

# (Andrew) Delete duplicated snps in 2 snp files
if ([ $JobType -eq 5 ]); then
  my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/remove_duplicated_snps.R
  R --vanilla --slave -f $my_script
fi

# Compute MAF per tissue
if ([ $JobType -eq 6 ]); then
  my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/compute_MAF_per_tissue.R
  R --vanilla --slave -f $my_script
fi

# Run MatrixEQTL by context
if ([ $JobType -eq 7 ]); then
  exp_scale=$2
  i=$3
  my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/run_MatrixEQTL_by_context.R
  R --vanilla --slave -f $my_script --args $exp_scale $i
fi


# Get number of SNPs tested per gene in each tissue
if ([ $JobType -eq 8 ]); then
  tissue_name=$2
  nSNPs=$3

  input_dir=/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/MatrixEQTL
  output_dir=/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/TreeQTL

  if ([ $nSNPs -eq 1 ]); then
    awk '{print $2}' ${input_dir}/${tissue_name}.v8.EUR.normalized_and_residualized_expression.all_pairs.txt | sort | uniq -c | awk '{print $2,$1}' > ${output_dir}/n_SNPs_per_gene_${tissue_name}.txt
  fi

  if ([ $nSNPs -eq 2 ]); then
  awk '{print $1}' ${input_dir}/${tissue_name}.v8.EUR.normalized_and_residualized_expression.all_pairs.txt | sort | uniq -c | awk '{print $2,$1}' > ${output_dir}/n_genes_per_SNPs_Heterogeneous.txt
  fi

  if ([ $nSNPs -eq 3 ]); then
    awk '{print $2}' ${input_dir}/${tissue_name}.v8.EUR.normalized_and_residualized_expression_homogeneous.all_pairs.txt | sort | uniq -c | awk '{print $2,$1}' > ${output_dir}/n_SNPs_per_gene_${tissue_name}.txt
  fi

  if ([ $nSNPs -eq 4 ]); then
  awk '{print $1}' ${input_dir}/${tissue_name}.v8.EUR.normalized_and_residualized_expression_homogeneous.all_pairs.txt | sort | uniq -c | awk '{print $2,$1}' > ${output_dir}/n_genes_per_SNPs_Homogeneous.txt
  fi


fi

# Filter MatrixEQTL results for eQTL p-value <= threshold
if ([ $JobType -eq 9 ]); then
  tissue_name=$2

  input_dir=/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/MatrixEQTL
  output_dir=/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/MatrixEQTL_FDRthreshold

  head -n 1 ${input_dir}/${tissue_name} > ${output_dir}/${tissue_name}
  awk '{if($5<=5e-01) print $0}' ${input_dir}/${tissue_name} >> ${output_dir}/${tissue_name}

fi

# Run TreeQTL
if ([ $JobType -eq 10 ]); then
  exp_scale=$2
  my_script=/u/project/zaitlenlab/bballiu/FastGxE/scripts/run_TreeQTL.R
  R --vanilla --slave -f $my_script --args $exp_scale
fi


# Summarize eAssociation files across tissues in GTEx
if ([ $JobType -eq 11 ]); then
    R --vanilla --slave -f /u/project/zaitlenlab/bballiu/FastGxE/scripts/summarize_treeqtl_output.R
fi

# Prepare input files for mashR (Andrew)
if ([ $JobType -eq 12 ]); then
    R --vanilla --slave -f /u/project/zaitlenlab/bballiu/FastGxE/scripts/create_mashR_input_files.R
fi

# Run mashR analysis
if ([ $JobType -eq 13 ]); then
    echo Under construction
fi

# eVariants enrichment in genomic features
if ([ $JobType -eq 14 ]); then

    R --vanilla --slave -f /u/project/zaitlenlab/bballiu/FastGxE/scripts/__main__run_tissue_agnostic_enrich_genomfeatures.R

fi
# Run Alexis Battle method (sn_spMF): make input files
if ([ $JobType -eq 15 ]); then

    R --vanilla --slave -f /u/project/zaitlenlab/bballiu/FastGxE/scripts/make_sn-spMF_input_files.R

fi

# Run Alexis Battle method (sn_spMF): model selection: grid search for hyperparameters (step 1: narrow down range)
if ([ $JobType -eq 16 ]); then

    echo "starting..."
    start=$(date +%s.%N)

    cd /u/project/zaitlenlab/bballiu/FastGxE/other_methods/sn-spMF/ts_eQTLs-master

    K=$2
    alpha1=$3
    lambda1=$4
    iterations=$5

    Rscript sn_spMF/run_MF.R \
    -k $K -a ${alpha1} -l ${lambda1} -t ${iterations} \
    -x "/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/sn-spMF/TbT_X.txt" \
    -w "/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/sn-spMF/TbT_SE.txt" \
    -O "/u/project/zaitlenlab/bballiu/FastGxE/results/eQTL_mapping/sn-spMF/Model_Selection/Step1_NarrowDownRange/"

    duration=$(echo "$(date +%s.%N) - $start" | bc)
    execution_time=`printf "%.2f seconds" $duration`
    echo "Script Execution Time: $execution_time"
    
    echo "ending..."
fi
