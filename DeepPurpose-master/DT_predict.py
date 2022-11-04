# --coding:utf-8--
# @Time : 2022/9/12 9:53
# @File : Predict.py
# @Software: PyCharm
import os
import re
import numpy as np
import pandas as pd
from DeepPurpose import DTI as models
from DeepPurpose import utils, dataset
from DeepPurpose.dataset import *





machine_target = pd.read_csv("../Network_Evaluation_Tools-master_2/result/Machine_target.csv",index_col=0)
target_seq = pd.read_csv("./data/human_seq_all.csv",index_col=0)
col_index = [r for r in list(machine_target) if re.search("_target",r)]
machine_target = machine_target[col_index]
machine_target = machine_target.dropna()
machine_target = machine_target.drop_duplicates()
machine_target = machine_target.loc[:,col_index]
drug = pd.read_csv("../input_smile/mxsgt_cf_smile.csv",index_col=0,encoding="GB2312")
drug = drug.iloc[:,0]
drug = drug.dropna()

t_index = []
for col in list(machine_target):
    m_t = machine_target[col]
    t_i = []
    for i in list(target_seq["Gene.names"]):
        b = i in list(m_t)
        t_i.append(b)
    t_index.append(t_i)

path = 'DrugBank_CNN_model'
net = models.model_pretrained(path_dir = path)
drug_encoding, target_encoding = 'CNN', 'CNN'

y_predict_result = pd.DataFrame()
for index in t_index:
    target_seq_1 = target_seq[index]
    #target_seq_1 = target_seq_1.drop_duplicates(subset="Gene.names",keep='first')
    target_seq_1.columns = ["target", "seq"]
    seq = target_seq_1['seq']
    tar = target_seq_1["target"]

    tar_1 = [tar]*len(drug)
    tar_1 = np.array(tar_1)
    tar_1 = tar_1.ravel(order = "F")

    t1 = [seq] * len(drug)
    t1 = np.array(t1)
    t1 = t1.ravel(order='F')
    d1 = list(drug) * len(seq)
    y = [-1] * len(d1)

    X_pred_1 = utils.data_process(d1, t1, y,
                                  drug_encoding, target_encoding,
                                  split_method='no_split')

    y_pred = net.predict(X_pred_1)
    y_preed_1 = pd.DataFrame(zip(y_pred, t1, d1,tar_1), columns=['pvalue', "target_sequence", "drug","target_names"])
    #y_screen = y_preed_1[y_preed_1['pvalue'] >= 0.5]
    y_screen = y_preed_1.sort_values(by="pvalue", ascending=False)
    y_predict_result = pd.concat([y_predict_result,y_screen],axis=1)

y_predict_result.to_csv("../result/y_predict_target.csv")

