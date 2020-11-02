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
library(SummarizedExperiment)
library(DESeq2)

cancers = c("ACC","BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC","ESCA",
            "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
            "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC",
            "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCS", "UVM")

query_project = function(cancer, sample_type = "Primary Tumor"){
  query_TCGA = GDCquery(
    project = paste0("TCGA-",cancer),
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "HTSeq - Counts",
    sample.type = c(sample_type))
  GDCdownload(query = query_TCGA)
  query = GDCprepare(query_TCGA)
  return(query)
}
reconcile_stage = function(s){
s = toupper(s)
  if(s %in% c("STAGE I", "STAGE IA","STAGE IB","STAGE IC")){
    return("Stage I")
  } else if(s %in% c("STAGE II", "STAGE IIA","STAGE IIB","STAGE IIC")){
    return("Stage II")
  } else if (s %in% c("STAGE III", "STAGE IIIA","STAGE IIIB","STAGE IIIC")){
    return("Stage III")
  }else if (s %in% c("STAGE IV", "STAGE IVA","STAGE IVB","STAGE IVC")){
    return("Stage IV")
  } else{
    return(s)
  }
}

normalize_counts = function(ht_counts, col_data){
  print("ping")
  deseq_counts=DESeqDataSetFromMatrix(ht_counts, colData = col_data,design = ~1)
  print("tang")
  deseq_counts = DESeq(deseq_counts)
  print("slang")
  normalized_counts = counts(deseq_counts, normalized=TRUE)
  normalized_counts = t(normalized_counts)
  normalized_counts = scale(normalized_counts)
  normalized_counts = t(normalized_counts)
  return(normalized_counts)
}


argv = commandArgs(trailingOnly = T)
cancer = argv[1]
sample_type = argv[2]

if(sample_type=="t"){
  sample_type = "Primary solid Tumor"
} else if(sample_type =="n"){
  sample_type = "Solid Tissue Normal"
}

if(cancer %in% cancers){ 
se = query_project(cancer=cancer,sample_type=sample_type)
dir.create(paste0("../data/input_data/tcga/",cancer),showWarnings = F)
dir.create(paste0("../data/input_data/tcga/",cancer,"/samples/named/"),recursive=T,showWarnings = F)
dir.create(paste0("../data/input_data/tcga/",cancer,"/samples/unnamed/"),recursive=T,showWarnings = F)
dir.create(paste0("../data/input_data/tcga/",cancer,"/samples/meta/"),recursive=T,showWarnings = F)

exp_data = assay(se)
patient_data = colData(se)

n_patients= ncol(exp_data)
set.seed(1234)
gene_meta = rowData(se)
for( i in seq(1,100)){
  subset = exp_data[,sample(n_patients, floor(0.9*n_patients),replace=F)]
 counts = normalize_counts(subset,colData(se)[colnames(subset),])
  fileConn<-file(paste0("../data/input_data/tcga/",cancer,"/samples/meta/",cancer,"_sample_",i,"names.txt"))
  writeLines(colnames(subset), fileConn)
  close(fileConn)
  write.table(data.frame(counts),
              paste0("../data/input_data/tcga/",cancer,"/samples/unnamed/",cancer,"_sample_",i,".txt"),
              row.names=F,col.names=F,sep=",")
   write.table(data.frame(counts),
              paste0("../data/input_data/tcga/",cancer,"/samples/named/",cancer,"_sample_",i,".txt"),
              row.names=T,col.names=F,sep=",")
}

counts = normalize_counts(exp_data,colData(se))
write.table(data.frame(counts),
              paste0("../data/input_data/tcga/",cancer,"/tumor_expression.txt"),
              row.names=F,col.names=F,sep=",")
write.table(data.frame(counts),
              paste0("../data/input_data/tcga/",cancer,"/named_tumor_expression.txt"),
              row.names=T,col.names=F,sep=",")
write.csv(data.frame(gene_meta),paste0("../data/input_data/tcga/",cancer,"/index_gene_map.csv"),row.names = T)

} else{
  cat("\nerror: cancer not a recognized project")
}
#unlink("GDCdata",recursive=TRUE)





