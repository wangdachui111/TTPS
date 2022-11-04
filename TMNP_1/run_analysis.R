rm(list = ls())

#检查是否安装了所需R包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")


Package_require <- c("dplyr","ggplot2","pheatmap","cowplot")



source("helper_TMNP.R")
source("help_plot_1.R")
source("result_return.R")





#读取输入数据
input_data <- dir("../input_rna_deg/")
input_data <- input_data[grep("RNAseq_DEG.csv",input_data)]
querydata <- read.csv(paste0("../input_rna_deg/",input_data),header = T)
if (is.character(querydata[,1])) {
  rownames(querydata) <- querydata[,1]
  querydata <- querydata[,-1]
}else if(is.numeric(querydata[,1])){
  querydata <- querydata
}


#
# result_dat <- fun_all(querydata = querydata)
# names(result_dat) <- c("Target_Gene","Pathway","GO_Biological_Process","Pathological_Process","Cell","Tissue","DisGenet")



if (Module == "Target Gene") {
  target <- dat_return(Module)
}else if (Module == "Pathway") {
  pathway <- dat_return(Module)
}else if (Module == "GO biological process") {
  go <- dat_return(Module)
}else if (Module == "Pathological process") {
  pp <- dat_return(Module)
}else if (Module == "Cell") {
  cell <- dat_return(Module)
}else if (Module == "Tissue") {
  tissue <- dat_return(Module)
}else if (Module == "DisGenet") {
  dis <- dat_return(Module)
}else if (Module == "all") {
  
}


Module = "Target Gene"
target <- dat_return(Module)
# target_plot <- target[[2]]
# target_plot

Module = "Pathway"
pathway <- dat_return(Module)
# pathway_plot <- pathway[[2]]
# pathway_plot

Module = "GO biological process"
go <- dat_return(Module)
# go_plot <- go[[2]]
# go_plot

Module = "Pathological process"
pp <- dat_return(Module)
# pp_plot <- pp[[2]]
# pp_plot

Module = "Cell"
cell <- dat_return(Module)
# cell_plot <- cell[[2]]
# cell_plot



Module = "Tissue"
tissue <- dat_return(Module)
# tissue_plot <- tissue[[2]]
# tissue_plot

Module = "DisGenet"
dis <- dat_return(Module)
# dis_plot <- dis[[2]]
# dis_plot

x = target[[1]]
plot_fun(res_all = target[[1]],Module = 1,res_name = colnames(querydata),p_name = "Target_Gene")





Module = "Tissue"
tissue <- tissue
plot_heatmap(res_all = tissue,Module = 1,res_name = colnames(querydata),p_name = "Tissue")
plot_hist(res_all = tissue,Module = 1,res_name = colnames(querydata),p_name = "Tissue")
plot_corr(res_all = tissue,Module = 1,res_name = colnames(querydata),p_name = "Tissue")

# # 清除图片缓存
# allfile = dir("../plot_result/")
# pdffile <- grep("*", allfile)
# file.remove(paste0("../plot_result/",allfile[pdffile]))

