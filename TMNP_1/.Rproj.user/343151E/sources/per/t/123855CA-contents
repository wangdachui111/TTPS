dat_return <- function(Module){
  if (Module == "Target Gene"){
    res_target <- fun_target(querydata)
    df <- res_target[[1]]
    write.csv(df,file = "../result/result_Target_Gene.csv")
    P <- plot_fun(res_all = res_target,Module = 1,res_name = colnames(querydata),p_name = "Target_Gene")
    return(list(res_target,P))
  }else if(Module == "Pathway"){
    res_pathway <- fun_Pathway(querydata)
    df <- res_pathway[[1]]
    write.csv(df,file = "../result/result_Pathway.csv")
    P <- plot_fun(res_all = res_pathway,Module = 1,res_name = colnames(querydata),p_name = "Pathway")
    return(list(res_pathway,P))
  }else if(Module == "GO biological process"){
    res_GO <- fun_GO(querydata)
    df <- res_GO[[1]]
    write.csv(df,file = "../result/result_GO_biological_process.csv")
    P <- plot_fun(res_all = res_GO,Module = 1,res_name = colnames(querydata),p_name = "GO_biological_process")
    return(list(res_GO,P))
  }else if(Module == "Pathological process"){
    res_pp <- fun_pp(querydata)
    df <- res_pp[[1]]
    write.csv(df,file = "../result/result_Pathological_process.csv")
    P <- plot_fun(res_all = res_pp,Module = 1,res_name = colnames(querydata),p_name = "Pathological_process")
    return(list(res_pp,P))
  }else if(Module == "Cell"){
    res_cell <- fun_cell(querydata)
    df <- res_cell[[1]]
    write.csv(df,file = "../result/result_Cell.csv")
    P <- plot_fun(res_all = res_cell,Module = 1,res_name = colnames(querydata),p_name = "Cell")
    return(list(res_cell,P))
  }else if(Module == "Tissue"){
    res_tissue <- fun_tissue(querydata)
    df <- res_tissue[[1]]
    write.csv(df,file = "../result/result_Tissue.csv")
    P <- plot_fun(res_all = res_tissue,Module = 1,res_name = colnames(querydata),p_name = "Tissue")
    return(list(df,P))
  }else if(Module == "DisGenet"){
    res_dis <- fun_disgenet(querydata)
    df <- res_dis[[1]]
    write.csv(df,file = "../result/result_DisGenet.csv")
    P <- plot_fun(res_all = res_dis,Module = 1,res_name = colnames(querydata),p_name = "DisGenet")
    return(list(res_dis,P))
  }
}





