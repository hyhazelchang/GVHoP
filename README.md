# GVHoP: Giant Virus Host Predictor

🖥️ GVHoP is a command-line tool for predicting eukaryotic host to giant virus metagenome-assembled genomes (GVMAGs).

GVHoP predicts hosts for GVMAGs by combining two complementary sources of genomic information. First, predicted viral proteins are screened against curated HMMs representing protein families. In parallel, viral proteins are compared with a eukaryotic protein database to identify potential virus–host gene similarities. These gene-content and gene-similarity signals are then integrated to infer the most likely eukaryotic hosts of the input GVMAGs.


## How to use

### Installation
Clone the `GVHoP` repository:

```{bash}
git clone https://github.com/hyhazelchang/GVHoP.git
```
Set up GVHoP

```{bash}
cd GVHoP/
chmod u+x gvhop
```

### Requirements
Create `conda` environment and install requirements:

```{bash}
conda env create -f env.yml
```

Download the databases from ZENODO before running GVHoP:

```{bash}
wget 'https://doi.org/10.5281/zenodo.22199562/GVHoP_database_v1.0.zip'
unzip GVHoP_database_v1.0.zip
```

### Example run

Activate GVHoP environment:

```{bash}
conda activate gvhop
```

Run GV-host predictor:

```{bash}
./gvhop --in_dir=example_MAGs/ --in_file_ext=faa --db_dir=GVHoP_database_v1.0/ --sample_ls=source_data/example_inputs/sample.ls --out_dir=gvhop_out/ --cpu=8
```

### Results

The output files includes:
=> Feature sets for prediction
   1.



#### Interpretation of results:

## Contact

