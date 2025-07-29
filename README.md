# aptamerSPECIAL
# aptamerSPECIAL项目功能与逻辑梳理讲义
## 项目概述
aptamerSPECIAL是一个用于生成、分析和分类适体（Aptamer）二级结构的完整计算分析 pipeline，主要通过模拟不同类型的适体序列、预测其二级结构并进行结构特征提取与机器学习分类，最终实现适体结构的自动化分析与可视化。

## 文件功能详解
### 1. 核心代码文件 `generate_aptamer_libraries.py`
- 功能 ：生成不同类型的适体序列库
- 关键函数 ：
  - generate_g_quadruplex() ：生成G-四链体结构序列（核心为G重复序列+环序列）
  - generate_pseudoknot() ：生成假结结构序列（茎环结构+交叉配对区域）
  - generate_triplex() ：生成三链结构序列（双链区域+Hoogsteen配对第三链）
  - generate_random_sequence() ：生成随机序列作为对照
  - main() ：创建输出目录，生成G-四链体、假结、三链体各2000条及随机序列3000条，保存至 `aptamer_libraries` 目录 `predict_secondary_structures.py`
- 功能 ：调用ViennaRNA Package的RNAfold工具预测RNA二级结构
- 关键流程 ：
  1. 1.
     解析命令行参数（输入FASTA文件、输出路径、RNAfold路径）
  2. 2.
     predict_secondary_structures() 函数：创建输出目录，构建并执行RNAfold命令
  3. 3.
     将预测结果保存至 `secondary_structures.txt` `analyze_aptamer_structures.py`
- 功能 ：结构特征提取、聚类分析与机器学习分类
- 核心模块 ：
  - 数据解析： parse_structure_file() 提取序列/结构/MFE/标签
  - 特征工程： calculate_structure_features() 计算长度/碱基对/发夹结构等特征
  - 聚类分析： cluster_sequences() 使用DBSCAN基于编辑距离聚类
  - 分类模型： train_classifier() 训练随机森林分类器并评估
  - 可视化：t-SNE/UMAP降维、特征分布箱线图、相关性热图等
### 2. 辅助文件 `requirements.txt`
- 功能 ：项目依赖管理
- 核心依赖 ：numpy、pandas、scikit-learn、python-Levenshtein（编辑距离计算）、umap-learn==0.5.3（降维） `run_analysis_pipeline.bat`
- 功能 ：Windows批处理脚本，自动化执行完整分析流程
- 执行步骤 ：
## 项目逻辑流程图
## 关键输出文件
1. 1.
   序列文件： `aptamer_libraries` （包含各类型FASTA文件）
2. 2.
   结构预测： `secondary_structures.txt`
3. 3.
   特征数据： `aptamer_features.csv` 、 `aptamer_with_predictions.csv`
4. 4.
   可视化结果：6种PNG格式图表（t-SNE/UMAP聚类图、特征热图等）
5. 5.
   分析报告： `analysis_report.txt`
## 使用说明
1. 1.
   安装依赖： pip install -r requirements.txt
2. 2.
   执行完整流程：双击 `run_analysis_pipeline.bat` 或命令行运行该脚本
3. 3.
   查看结果：所有输出文件将保存在项目根目录下
