#DAISY Lineage Tracing Analysis
This is the preprocessing step of amplicon sequence sublibrary for lineage tracing.  
We use CRISPR/Cas12a-based lineage tracing named DAISY published in Mol Cell 2022 by Dr. Le Cong lab.  
Please change your path and file for your NGS output fastq file.  


###2024/7/9 Upload LineageTracing_MO1


###2025/2/18 Upload LineageTracing_MO2 for cassiopeia input


##Required envinronment.  
Could you set conda environment as follows before use, please?

###For DAISY preprocessing,

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

###For Cassiopeia analysis,

conda create -n cassiopeia python==3.10.13
conda activate cassiopeia

pip install git+https://github.com/YosefLab/Cassiopeia@master#egg=cassiopeia-lineage

