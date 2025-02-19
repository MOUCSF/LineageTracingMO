# DAISY Lineage Tracing Analysis
This is the preprocessing step of amplicon sequence sublibrary for lineage tracing.  
Our pipeline can analyze staticBC (intBC) and mutableBC from .fastq data to .txt which describes clonal population with mutation information.
We use CRISPR/Cas12a-based lineage tracing named **DAISY** published in Mol Cell 2022 by Dr. Le Cong lab.  
We also utilize **Cassiopeia** published in Genome Biol 2020 and **LAML** published in Biorxiv 2024 for analysis.
Please cite their papers if you use for your works.  
Nicholas W. Hughes, Yuanhao Qu, Jiaqi Zhang. et al. Machine-learning-optimized Cas12a barcoding enables the recovery of single-cell lineages and transcriptional profiles, Molecular Cell, 82, 16 (2022). https://doi.org/10.1016/j.molcel.2022.06.001  
Jones, M.G., Khodaverdian, A., Quinn, J.J. et al. Inference of single-cell phylogenies from lineage tracing data using Cassiopeia. Genome Biol 21, 92 (2020). https://doi.org/10.1186/s13059-020-02000-8  
Uyen Mai, Gillian Chu, Benjamin J. Raphael. et al. Maximum Likelihood Inference of Time-scaled Cell Lineage Trees with Mixed-type Missing Data. BioRxiv (2024). https://doi.org/10.1101/2024.03.05.583638  


### 2024/7/9 Upload LineageTracing_MO1 for iqtree2 input
I used all genome editing including CRISPR mediated mutagenesis as well as indel for generating 120 character (120 bp) matrix.

### 2025/2/18 Upload LineageTracing_MO2 for cassiopeia input
I used only Indel information for cassiopeia analysis for generating 4 character (CRISPR target 30 bp is 1 character after CIGAR analysis) x copy number. 

## Required envinronment.  
Could you set conda environment as follows before use, please? I am not familiar with dependency.

### For DAISY preprocessing,
```
conda create -n daisy python==3.10.13
conda activate daisy

conda install -c bioconda fastp
conda install -c conda-forge rapidfuzz 
conda install -c conda-forge scikit-bio
conda install conda-forge::biopython
conda install bioconda::emboss

python3 -m pip install -U dendropy
conda install -c bioconda iqtree
pip install laml
pip install edlib
```
### For Cassiopeia analysis,
```
conda create -n cassiopeia python==3.10.13
conda activate cassiopeia

pip install git+https://github.com/YosefLab/Cassiopeia@master#egg=cassiopeia-lineage
```
## How to use.
This tutorial is mainly for ver2 code. Please change path of your working folder and files for your NGS output fastq file after downloading my codes.  

Before use, 10x barcode whitelist files are required.  
For single cell barcode which I used as whilelist in the code, first please make seurat object by using your cellranger output as usual.  
After removal of low quality cells, dead cells, and doublet cells by your familiar fashion, you can extract 10x barcode and save as yoursample.whitelist.Singlet.txt as follows. My case, I used SampleID as obj@meta.data$SampleID to distinguish sample number.  
Please type in **Rstudio** as follows if you use merged seurat object.
```
save_10xbarcode <- function(obj) {
  SampleIDs <- unique(obj$SampleID) # SampleID の一覧を取得
  for (SampleID in SampleIDs) {
    cells <- colnames(obj)[obj$SampleID == SampleID] # 該当するSampleIDのセルバーコードを取得
    tenxbarcode <- sub("-1$", "", cells) # "-1" を削除
    tenxbarcode <- gsub(paste0("^", SampleID, "_"), "", tenxbarcode) # "SampleID_" の部分を空白に置き換えて削除
    write.table(tenxbarcode, file = paste0(SampleID, ".whitelist.singlet.txt"),quote = FALSE, row.names = FALSE, col.names = FALSE)}# ファイルに保存
  }
save_10xbarcode(obj)
```
Please type in **Rstudio** as follows if you have each single seurat object.
```
save_10xbarcode <- function(sample_name){
  obj <- get(sample_name, envir = .GlobalEnv)  # sample_nameに対応するオブジェクトを取得
  tenxbarcode_with1 <- colnames(obj)
  tenxbarcode <- sub("-1$", "", tenxbarcode_with1) # "-1" を削除
  write.table(tenxbarcode, file=paste0(sample_name, ".whitelist.Singlet.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
}

sample_name_list <- c("DAISY0","DAISY13","DAISYRETRY") # Your seurat object name.
lapply(sample_name_list, save_10xbarcode)
```
After you set yoursample.whitelist.Singlet.txt in appropriate path, please copy and paste command as follows in the **terminal** or **linux** in your daisy environment.
```
bash /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step0ver1.sh && bash /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step1ver1.sh && bash /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step2ver3.sh && python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step3ver2.py && python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step4.py && python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step5.py && python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step7ver3.py && python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step8ver5.py && python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step9.py
```
Depending on your sublibrary sequencing file size, it takes one week to complete entire analysis of fastq.gz data. Please be patient.  
And then, please run Cassiopeia in your cassiopeia environment. This step is not set for all samples. Please change sample number one by one in your code.
```
python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step11.py
```
or
```
python /Users/masahirookada/Desktop/Guptadaisy/mo12398/mo12398_step12.py
```

We will update them for more user friendly manner. I hope you can run successfully.


Sincerely,  
Masahiro Okada  
This work is supported by Dr. Diaz and Dr. Gupta.
