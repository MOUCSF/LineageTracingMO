# DAISY Lineage Tracing Analysis
This is the preprocessing step of amplicon sequence sublibrary for lineage tracing.  
We use CRISPR/Cas12a-based lineage tracing named DAISY published in Mol Cell 2022 by Dr. Le Cong lab.  
We also utilize Cassiopeia published in Genome Biol 2020 and LAML published in Biorxiv 2024.
Please cite their papers if you use.  
Jones, M.G., Khodaverdian, A., Quinn, J.J. et al. Inference of single-cell phylogenies from lineage tracing data using Cassiopeia. Genome Biol 21, 92 (2020). https://doi.org/10.1186/s13059-020-02000-8  
Nicholas W. Hughes, Yuanhao Qu, Jiaqi Zhang. et al. Machine-learning-optimized Cas12a barcoding enables the recovery of single-cell lineages and transcriptional profiles, Molecular Cell, 82, 16 (2022). https://doi.org/10.1016/j.molcel.2022.06.001  
Uyen Mai, Gillian Chu, Benjamin J. Raphael. et al. Maximum Likelihood Inference of Time-scaled Cell Lineage Trees with Mixed-type Missing Data. BioRxiv (2024). https://doi.org/10.1101/2024.03.05.583638  

### 2024/7/9 Upload LineageTracing_MO1


### 2025/2/18 Upload LineageTracing_MO2 for cassiopeia input


## Required envinronment.  
Could you set conda environment as follows before use, please?

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
Please change your path and file for your NGS output fastq file.  



