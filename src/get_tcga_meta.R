##
#
#	 process the metadata of tcga to help choose experimental cancers
#				
#
###

library(TCGAbiolinks)
cancers = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP",
	"LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD",
	"TGCT","THCA","THYM","UCS","UVM")

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



	
for (cancer in cancers){
	pjct = paste0("TCGA-",cancer)
	fpath= paste0("../data/input_data/tcga/",cancer)
	dir.create(fpath)
	query.exp <- GDCquery(project =pjct , 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))
	GDCdownload(query.exp)

	#expr_data <- GDCprepare(query = query.exp, save = TRUE, save.filename = "brcaExp.rda")
	# Which samples are Primary Tumor
	dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP")
	# which samples are solid tissue normal
	dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")
	dataClin <- GDCquery_clinic(project = pjct,"clinical") 

	if("ajcc_pathologic_stage" %in% colnames(dataClin)){
		dataClin["Stage"]=apply(dataClin["ajcc_pathologic_stage"],1,reconcile_stage)
	} else {
		dataClin["Stage"]=apply(dataClin["tumor_stage"],1,reconcile_stage)
	}
	
	
	n_tumor_samples = length(dataSmTP) 
	n_normal_samples = length(dataSmNT)
	stage_table = table(dataClin["Stage"])
	dir.create(paste0("../data/input_data/tcga/",cancer,"/meta/"),recursive=T,showWarnings=F)
	fileConn<-file(paste0("../data/input_data/tcga/",cancer,"/meta/n_samples.txt"))
  	writeLines(paste0('Tumor:\t',n_tumor_samples ,'\n','Normal:\t',n_normal_samples), fileConn)
  	write.csv(stage_table,paste0("../data/input_data/tcga/",cancer,"/stage_counts.csv"),row.names = T)
  	close(fileConn)  
}
	