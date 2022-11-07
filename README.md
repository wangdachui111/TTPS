# TTPS
Transcriptome-based target prediction system
Required Data
You can download the data required for caculation on 
## Installation
### TTPS相关运行程序需要用到三种运行环境：*TMNP_1(R4.2.0)*，*Network_Evaluation_Tools-master_2(python2.7)*，*DeepPurpose-master(python3.6)*。
#### TMNP_1 
需要安装4.2.0以上版本，所需R包可以运行*TMNP_1/packages_installation.R*
#### Network_Evaluation_Tools-master_2
```
conda create -n py27 python=2.7
conda activate py27
pip install argparse
pip install neworkx
pip install scipy
pip install scikit-learn
```
#### DeepPurpose-master
```
conda create -n DeepPurpose python=3.6
conda activate DeepPurpose
conda install -c conda-forge rdkit
pip install descriptastorus 
pip install DeepPurpose
```


## 简介
TTPS是一个具有多种生物信息分析手段的靶点预测系统，致力于预测药物诱导的不同层次的生物效应。该方法集成了多个高通量数据库的信息，我们将药物的生物效应在靶点、通路、生物过程、病理过程、细胞水平、组织水平、疾病反应等七个生物水平进行系统划分，我们的预测系统可以将药物在不同系统的干扰效应通过NCS值打分来量化，并且应用了统计学方法对生物效应的显著性进行评估。转录组数据可以揭示药物在全基因水平上的变化，但是大量的基因靶点之间彼此关联，我们的预测系统对于预测靶点进一步地区分，可以预测药物的直接靶点和间接靶点。我们将两种不同类型的靶点映射到通路互作网络中，通过社区模块化方法，我们可以借此发现靶点在通路模块的关联。靶点由直接靶点到间接靶点的作用方向可以帮助我们更加清楚地追踪靶点的效应，结合通路的分子机制可以帮助我们辨别复杂的通路效应。
该工具使用了R语言(R4.2.0)和python语言(python2.7;python3.6)。提供了如下子工具：
* 药物与--个疾病反应的NCS打分和显著性评估(Module = "Target Gene")
* 药物在32个组织上的影响效应的NCS打分和显著性评估(Module = "Target Gene")
* 药物在3021种细胞类型上的影响效应的NCS打分和显著性评估(Module = "GO biological process")
* 药物在－－种病理过程中的影响效应的NCS打分和显著性评估(Module = "Pathological process")
* 药物在－－个生物过程上的影响效应的NCS打分和显著性评估(Module = "Cell")
* 药物在－－个通路水平上的影响效应的NCS打分和显著性评估(Module = "Tissue")
* 药物在3275个基因靶点上的影响效应的NCS打分和显著性评估(Module = "DisGenet")
* 通过药物的smile和靶点的氨基酸结构预测药物和靶点的直接结合概率


## 用法说明
对于预测药物对七个生物水平的影响效应，我们应用R程序进行了预测，并简单绘制了相应的数据分布图。运行*TMNP_1/run_analysis.R*文件，设置不同的*Module*参数,即可计算相应的模块打分和显著性评估，相应的数据和图分别保存在*result*和*plot_result*文件夹下。对于药物和基因靶点直接结合概率的预测，我们首先需要运行*Network_Evaluation_Tools-master_2/Machine_target.py*文件，随后运行*DeepPurpose-master/DT_predict.py*文件，结果保存为*result/y_predict_target.csv*。
