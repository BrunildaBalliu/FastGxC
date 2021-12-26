#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%% Script to decompose expression into shared and context specific components
#%%%%%%%%%%%%%%% Brunilda Balliu 
#%%%%%%%%%%%%%%% December 1st 2021
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%  Directories for input and output files 
args=commandArgs(TRUE)
work_dir = args[1]
exp_mat_filename = args[2]
data_dir=paste0(work_dir,'data/') # directory with phenotype and genotype data1
if(!dir.exists(data_dir)) dir.create(data_dir)
source(paste0(work_dir, "functions.R"))

#%%%%%%%%%%%%%%% R libraries 
library(data.table)
library(reshape2)
library(magrittr)

#%%%%%%%%%%%%%%% Read expression matrix, genes in columns, samples in rows.
exp_mat=read.table(file = paste0(data_dir,exp_mat_filename), sep = '\t')

#%%%%%%%%%%%%%%% Print number of genes and samples
sprintf("There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(exp_mat), ncol(exp_mat),max(colSums(is.na(exp_mat))),max(rowSums(is.na(exp_mat))))

#%%%%%%%%%%%%%%% Sample and context names
design=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(rownames(exp_mat)[i], split = "_"))[1])
context_names=sapply(1:nrow(exp_mat), function(i) unlist(strsplit(rownames(exp_mat)[i], split = "_"))[2])
contexts=unique(context_names)

#%%%%%%%%%%%%%%% Decompose expression into homogeneous and heterogeneous context expression
print("Decomposing data")
dec_exp_all=decompose(X = exp_mat, design = design)
bexp_all=dec_exp_all$Xb
wexp_all=dec_exp_all$Xw
bexp_all[is.nan(bexp_all)]=NA
wexp_all[is.nan(wexp_all)]=NA

sprintf("Between individual matrix: There are %s individuals and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(bexp_all), ncol(bexp_all),max(colSums(is.na(bexp_all))),max(rowSums(is.na(bexp_all))))

sprintf("Within individual matrix: There are %s samples and %s genes. The max number of missing samples for a gene is  %s. The max number of missing genes for a sample is  %s.", nrow(wexp_all), ncol(wexp_all),max(colSums(is.na(wexp_all))),max(rowSums(is.na(wexp_all))))

#%%%%%%%%%%%%%%% Save decomposed expression files 
print("Finished decomposition, saving files")

print("Saving between-individuals variation matrix")
fwrite(x = data.table(t(bexp_all),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}, 
       file = paste0(data_dir,"context_shared_expression.txt"), quote = F, row.names = F, 
       col.names = T, append = F, sep = '\t')

print("Saving within-individuals variation matrix for context: ")
for(i in 1:length(contexts)){
  print(contexts[i])
  wexp_t = wexp_all[grep(pattern = paste0(contexts[i],"$"), rownames(wexp_all)),]
  rownames(wexp_t)=gsub(pattern = paste0("_",contexts[i]), replacement = "", x = rownames(wexp_t))
  fwrite(x = data.table(t(wexp_t),keep.rownames = T) %>% {setnames(., old = "rn", new = "geneID")[]}, 
         file = paste0(data_dir,contexts[i],"_specific_expression.txt"),quote = F, row.names = F, 
         col.names = T, append = F, sep = '\t')
}

