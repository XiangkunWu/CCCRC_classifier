#load R package-----------------------
options(stringsAsFactors  = F)
suppressMessages(library(tidyverse))
suppressMessages(library(getopt))
suppressMessages(library(caret))
suppressMessages(library(xgboost))
suppressMessages(library(tibble))

#set the working path------------------
#setwd("")
#load parameters--------------------
#在这里，加载对应的参数：
#Note: The input data of the model needs to be Log2 converted and scale standardized, and the user can judge and modify it according to the actual format of the data
#RNA-seq platform data with TPM or FPKM should be log2(expr+1) and scale standardization.
#RNA-seq platform data with voom should be scale standardization.
#AFFY platform data with RMA should be scale standardization.
args <- list(input.data= "CCCRC_CRCAFFYcohort_GEP.txt",   #expression profile file
          output.data = "CCCRC_output.txt", #output file
          log2 = "F",        # Whether the data needs to be Log2 processed:T/F,Default is T
          scale = "T",      # character","Whether the data needs to be scaled:T/F,Default is T
          xgboost.model = "CCCRC_classifier-model.Rdata", # model file  ：CCCRC_classifier-model.Rdata
          genelist = "CCCRC_classifier-genelist.Rdata")    # gene list：CCCRC_classifier-genelist.Rdata
#Load the predictive model and  gene list ---------------
xgboost_model <- readRDS(args$xgboost.model) #the CCCRC classifier
Genelist <- readRDS(args$genelist) #80 model genes

#Load data------------
#Note: modify the user's own input data here, the data format is gene x sample, 
#the first column is set to "Gene",please refer to "test_data.txt" for the detailed format; 
predict_data <- read.table(args$input.data,sep = "\t",header = T,check.names = F)
predict_data <- filter(predict_data,Gene %in% Genelist)
predict_data <- column_to_rownames(predict_data,"Gene")

#Check if the data covers the genes required by the model------------
#Note: If there are more than 2 missing genes, an error will be returned; if the missing genes are within 2 or less, it will be filled by the knn algorithm
miss_gene <- setdiff(Genelist,rownames(predict_data))
if (length(miss_gene) > 2) {
  cat(paste0("Missing Gene:",paste0(miss_gene,collapse = ",")))
  q(status = 1)
} else if (length(miss_gene) >0 & length(miss_gene) <=2 ){
  library(DMwR) #remotes::install_github("cran/DMwR")
  predict_data[miss_gene,] <- NA
  predict_data <- knnImputation(predict_data) 
} 
#Data standardization: perform log2(expr+1) and scale standardization---------
#Note: The input data of the model needs to be Log2 converted and scale standardized, and the user can judge and modify it according to the actual format of the data
#RNA-seq platform data with TPM or FPKM should be log2(expr+1) and scale standardization.
#RNA-seq platform data with voom should be scale standardization.
#AFFY platform data with RMA should be scale standardization.
is.log2 <- ifelse(is.null(args$log2),T,as.logical(args$log2))
is.scale <- ifelse(is.null(args$scale),T,as.logical(args$scale))
if (is.log2) {
  predict_data <- log2(predict_data+1)
}
if (is.scale) {
  predict_data <- t(scale(t(predict_data))) %>% data.frame(stringsAsFactors = F)
}
#prediction-------------
predict_res <- data.frame(id = colnames(predict_data),
                          CCCRC_predict = predict(object = xgboost_model,newdata = t(predict_data),type = "raw"))
write.table(predict_res,args$output.data,sep = "\t",quote = F,row.names = F)
#-------------
cat("Finish!")
q(status = 0)


###model evaluate
library(caret)
data1<-read.table("CCCRC_CRCAFFYcohort.txt",header = T,sep = "\t",check.names = F)
data2<-read.table("CCCRC_output.txt",header = T,sep = "\t",check.names = F)
obs<-as.factor(data1$CCCRC)
prediction<-as.factor(predict_res$CCCRC_predict)
confusionMatrix(prediction,obs)

