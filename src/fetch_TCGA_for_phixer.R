####
#
#   Code to fetch data from GDC/TCGA projects and
#     preprocess it for use in the phixer algorithm 
# 
#               @ James Bannon
#####


## Basic Pipeline: 
#
#   Input: A cancer project shorthand from the command line (see cancers vector)
#   Output: Five files, a meta data with counts for each cancer and
#           the following pair for both tumor and normal:
#               - m genes x n samples
#               - mapping index to various gene symbols (e.g. rowData) 
# 
#   Process:
#           -> read in cancer project abbreviation
#           -> run GDC query and download data
#           -> subset to just the samples with matched normal and tumor
#           -> split into different data sets
#           -> normalize the counts of the separate data sets with DESeq2
#           -> center and scale each data set
#           -> write the above files

library(TCGAbiolinks)
library(edgeR)
library(SummarizedExperiment)
library(DESeq2)
library(rDGIdb)
cancers = c("ACC","BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC","ESCA",
            "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
            "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC",
            "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCS", "UVM")

query_project = function(cancer){
  query_TCGA = GDCquery(
    project = paste0("TCGA-",cancer),
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "HTSeq - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal"))
  GDCdownload(query = query_TCGA)
  query = GDCprepare(query_TCGA)
  return(query)
}

stratify_samples = function(se,cancer){
  colinfo = colData(se)
  n_patients = length(unique(colinfo$patient))
  matched_patients = colinfo$patient[which(duplicated(colinfo$patient))]
  subset = colinfo[colinfo$patient %in% matched_patients,]
  n_matched = length(unique(subset$patient))
  write(paste0("num patients:\t ",n_patients,"\nnum matched:\t ",n_matched),
        paste0("../data/input_data/tcga/",cancer,"/meta.txt"))
  normal_samples = rownames(subset[subset$sample_type=="Solid Tissue Normal",])
  tumor_samples = rownames(subset[subset$sample_type=="Primary Tumor",])
  
  normal_data = colinfo[normal_samples,]
  tumor_data = colinfo[tumor_samples,]
  return(list(ns=normal_samples,nd=normal_data,
              ts=tumor_samples,td=tumor_data))

}

reconcile_genes = function(se, normal_samples, 
                           tumor_samples,cancer){
  exp_data = assay(se) 
  normal_expr = exp_data[,normal_samples]
  tumor_expr = exp_data[,tumor_samples]
  keep_normal = filterByExpr(normal_expr,min.count=10)
  keep_tumor = filterByExpr(tumor_expr,min.count=10)
  genes_in_both = intersect(names(keep_normal), names(keep_tumor))
  write(paste0("number of kept genes:\t ",length(genes_in_both)),
        paste0("../data/input_data/tcga/",cancer,"/meta.txt"),append = T)
  normal_expr = normal_expr[genes_in_both,]
  tumor_expr = tumor_expr[genes_in_both,]
  return(list(ne=normal_expr, te=tumor_expr,kept=genes_in_both))
}

normalize_counts = function(ht_counts, col_data){
  deseq_counts=DESeqDataSetFromMatrix(ht_counts, colData = col_data,design = ~1)
  deseq_counts = DESeq(deseq_counts)
  normalized_counts = counts(deseq_counts, normalized=TRUE)
  normalized_counts = t(normalized_counts)
  normalized_counts = scale(normalized_counts)
  normalized_counts = t(normalized_counts)
  return(normalized_counts)
}

write_files = function(tcounts, ncounts, meta,cancer){
  write.table(data.frame(tcounts),
              paste0("../data/input_data/tcga/",cancer,"/tumor_expression.csv"),
              row.names=F,col.names=F,sep=",")
  num_tumor_samples = dim(tcounts)[2]
  num_tumor_genes = dim(tcounts)[1]
  
  fileConn<-file(paste0("../data/input_data/tcga/",cancer,"/tdims.txt"))
  writeLines(paste0(num_tumor_samples,'\n',num_tumor_genes), fileConn)
  close(fileConn)  
  write.table(data.frame(ncounts),
              paste0("../data/input_data/tcga/",cancer,"/normal_expression.csv"),
              row.names=F,col.names=F,sep=",")
  
  num_normal_samples = dim(ncounts)[2]
  num_normal_genes = dim(ncounts)[1]

  fileConn<-file(paste0("../data/input_data/tcga/",cancer,"/ndims.txt"))
  writeLines(paste0(num_normal_samples,'\n',num_normal_genes), fileConn)
  close(fileConn) 

  write.csv(meta,paste0("../data/input_data/tcga/",cancer,"/index_gene_map.csv"),row.names = T)
}
argv = commandArgs(trailingOnly = T)
cancer = argv[1]
if(cancer %in% cancers){ 
query = query_project(cancer=cancer)
print("query done")
dir.create(paste0("../data/input_data/tcga/",cancer),showWarnings = F)
samples = stratify_samples(query,cancer=cancer)
reconciled_samples= reconcile_genes(query,normal_samples = samples$ns,
                                    tumor_samples = samples$ts,
                                    cancer=cancer)
tumor_counts = normalize_counts(reconciled_samples$te,col_data = samples$td)
normal_counts = normalize_counts(reconciled_samples$ne,col_data = samples$nd)
gene_meta = rowData(query)[reconciled_samples$kept,]
gene_meta = data.frame(gene_meta)
rownames(gene_meta) = c()
write_files(tumor_counts,normal_counts,gene_meta,cancer)
} else{
  cat("\nerror: cancer not a recognized project")
}





