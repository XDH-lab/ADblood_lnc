#  Association of blood-based DNA methylation of lncRNAs with Alzheimer's disease diagnosis**

Guanghui Tian, Dehua Zheng, Nihui Zhang, Wanxin Deng, Jiankai Xu, Chenqing Huang, Qiang Chen, Shasha Hu, Li Xin, Hong Wang, Bo Wang, Kongning Li, Dahua Xu


# Description

Logistic regression analysis was conducted on DNA methylation data derived from blood samples of patients with Alzheimer’s disease (AD) and cognitively normal controls to identify epigenetically regulated (ER) long non-coding RNAs (lncRNAs). Using five machine learning algorithms, ER lncRNAs associated with AD diagnosis were prioritized. A blood-based diagnostic model for AD was subsequently developed based on lncRNA methylation profiles in participants from the AIBL study and independently validated in two large-scale, blood-derived cohorts: the ADNI and the European AddNeuroMed consortium.

### 1.Logistic Regression


|File            |Dataset                        |Link                         |
|----------------|-------------------------------|-----------------------------|
|Logistic Regression.R|`AIBL&ADNI&AddNeuroMed`            |      https://github.com/XDH-lab/ADblood_lnc/blob/main/Logistic%20Regression/Logistic%20regression.R    


### 2.Machine learning
|File            |Dataset                        |Link                         |
|----------------|-------------------------------|-----------------------------|
|Boruta.R        |`AIBL`            |https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/Boruta.R|
|LassoLR.R       |`AIBL`            |https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/LassoLR.R|
|Pamr.R          |`ABIL`|https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/Pamr.R|
|rfFuncs.R       |`AIBL`            |https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/rfFuncs.R|
|Xgboost.R       |`ABIL`|https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/Xgboost.R|
|venn.R          |`-`|https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/venn.R|
|Logistic Regression Modeling.R|`ABIL&ADNI&AddNeuroMed`|https://github.com/XDH-lab/ADblood_lnc/blob/main/ML/Logistic%20Regression%20Modeling.R|

### 3.Nomogram
|File            |Dataset                        |Link                         |
|----------------|-------------------------------|-----------------------------|
|Nomogram.R      |`AIBL&ADNI&AddNeuroMed`            |https://github.com/XDH-lab/ADblood_lnc/blob/main/Nomogram/Nomogram.R    





# For reproducible research
The following R packages are required:
```r
library(dplyr)
library(tidyverse)
library(minfi)
library(wateRmelon)
library(SummarizedExperiment)
library(ExperimentHub)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(caret)
library(randomForest)
library(glmnet)
library(xgboost)
library(pamr)
library(Boruta)
library(pROC)
library(rms)
library(car)
library(ggplot2)
library(venn)
library(VennDiagram)
library(doParallel)
library(foreach)
```
For ADNIMERGE, download it from [https://ida.loni.usc.edu/](https://ida.loni.usc.edu/): Merged ADNI 1/GO/2 Packages for R
```r
install.packages("/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
```
The platform information are:
```r
version   R version 4.3.2 (2021-05-18)
 os       CentOS Linux 8          
 system   aarch64, linux4.18​          
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
```
# Acknowledgement
Data used in preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in analysis or writing of this report. A complete listing of ADNI investigators can be found at: [http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf](http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf)
