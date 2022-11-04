# --coding:utf-8--
# @Time : 2022/8/19 16:47
# @File : Sys evaluation.py
# @Software: PyCharm
import os
import re

import pandas as pd
import networkx as nx
import numpy as np
from network_evaluation_tools import data_import_tools as dit
from network_evaluation_tools import network_evaluation_functions as nef
from network_evaluation_tools import network_propagation as prop

#读取过滤后的分子网络和基因集
network = dit.load_network_file('./Data/PID_Symbol.sif', verbose=True)
result_dir = os.listdir("../result/")
path_dir = [r for r in result_dir if re.search("target",r)]
df = pd.read_csv(path_dir,index_col=0)
df = df.iloc[:,1:]
group = df.shape[1]/4
node_fdr_all = df.iloc[:, group*3:]
nodesets_name = list(node_fdr_all)
node_fdr_sets = []
for node_name in nodesets_name:
    node_fdr = node_fdr_all[node_name]
    node_fdr_1 = node_fdr[node_fdr <= 0.05]
    node_fdr_set = [node_name]
    node_fdr_set.extend(node_fdr_1.index)
    node_fdr_sets.append(node_fdr_set)
genesets = {node[0]:set(node[1:]) for node in node_fdr_sets}

node_ncs_value = df.iloc[:,group:group*2]
node_ncs_value.columns = nodesets_name

#计算基因集子样本率
genesets_p = nef.calculate_p(network, genesets)
#确定最优传播常数
#alpha = prop.calculate_alpha(network)
alpha = 0.587
#计算传播的网络核心
kernel = nef.construct_prop_kernel(network, alpha=alpha, verbose=True)
#计算每个基因组的AUPRC值

#返回预测的扩展靶点
F_score = nef.small_network_Fscore_wrapper(kernel, genesets, genesets_p, n=30, cores=1, verbose=True)
F_score.columns = [geneset for geneset in genesets]

Machine_target = pd.DataFrame()
key_index = list(F_score)
key_value = [F_score[key].dropna() for key in F_score]
for i in range(len(key_index)):
    key_ncs = []
    key_value_1 = list(key_value[i])
    col_index = key_index[i]
    for node in key_value_1:
        node_neighbor = [n for n in network.neighbors(node)]
        node_value = node_ncs_value.loc[node_neighbor, col_index]
        node_value = node_value.dropna()
        node_value = [abs(n) for n in node_value]
        node_value_ncs = np.mean(node_value)
        key_ncs.append(node_value_ncs)
    key_value_1.extend(list(node_ncs_value.loc[genesets[col_index], col_index].index))
    key_ncs.extend(list(node_ncs_value.loc[genesets[col_index], col_index]))
    node_key_name = pd.Series(key_value_1)#[key_value_1, set(node_ncs_value.loc[key_value_1, key_index[i]].index)]
    node_key_ncs = pd.Series(key_ncs)#[key_ncs, set(node_ncs_value.loc[key_value_1, key_index[i]])]
    target = pd.concat([node_key_name, node_key_ncs], axis=1)
    target.columns = [col_index[:-3] + "target",col_index[:-3] + "NCS"]
    Machine_target = pd.concat([Machine_target, target],axis=1)

Machine_target.to_csv("Machine_target.csv")
