# GVHoP: Giant Virus Host Predictor

🖥️ GVHoP is a command-line tool for predicting eukaryotic host to giant virus metagenome-assembled genomes (GVMAGs) using gene content and giant virus-eukaryotes gene similarities

## How to use

### Installation
Clone the `GVHoP` repository:

```{bash}
git clone https://github.com/hyhazelchang/GVHoP.git
```

### Requirements
Create `conda` environment and install requirements:

```{bash}
conda env create -c conda-forge -c bioconda --name gvhop --file env.yml
```

Download the databases from ZENODO before running GVHoP:

```{bash}
wget '10.5281/zenodo.22199562' -O GVHoP_database_v1.0.zip
unzip GVHoP_database_v1.0.zip
```

### Example run

#### Activate GVHoP environment:

```{bash}
conda activate gvhop
```

#### Construct feature set (Gene content) from sequence data:

📝 This will generate the shell scripts for hmmsearch execution in sh_dir
```{bash}
python ExtractFeatureSets/cmd_hmm.py --in_dir=example_MAGs/ --db_dir=GVHoP_database_v1.0/GVHoP_GVOGs.hmm --out_dir=temp_out/output_hmm/ --sh_dir=temp_out/sh/hmmsearch/ --in_file_ext=faa --out_file_ext=domout --opt='hmmsearch --cpu 10 -E 1e-5 --domtblout' --n_job=4
```

🏃 Execute hmmsearch
```{bash}
python ExtractFeatureSets/execute_sh.py --sh_dir=temp_out/sh/hmmsearch/ --log_dir=temp_out/execute/hmmsearch/
```

💻 Parse feature matrix
```{bash}
python ExtractFeatureSets/parse_hmmsearch.py --in_dir=temp_out/output_hmm/ --in_file_ext=domout --out_dir=source_data/example_inputs/ --out_preffix=ex_GVOGs --cogset=source_data/features/GVHoP_GVOGs_all.tsv --seq_dir=example_MAGs/ --seq_file_ext=faa
```

#### Construct feature set (GV-euk signals) from sequence data:

📝 This will generate the shell scripts for DIAMOND blastp execution in sh_dir
```{bash}
python3 ExtractFeatureSets/cmd_diamond.py --task=blastp --in_dir=example_MAGs/ --out_dir=temp_out/output_blastp/ --sh_dir=temp_out/sh/blastp/ --in_file_ext=faa --out_file_ext=txt --db_dir=GVHoP_database_v1.0/GVHoP_GVEUKs.dmnd --outfmt=6 --threshold=8 --opt="--evalue 1e-5" --n_job=8
```

🏃 Execute diamond blastp
```{bash}
python ExtractFeatureSets/execute_sh.py --sh_dir=temp_out/sh/blastp/ --log_dir=temp_out/execute/blastp/
```

💻 Parse feature matrix
```{bash}
python ExtractFeatureSets/parse_blasthits.py --blast_in=temp_out/output_blastp/ --in_file_ext=txt --column_names=source_data/features/GVHoP_GVEUKs_all.tsv --outfile=source_data/example_inputs/ex_GVEUKs.tsv
```

Output feature sets for GVHoP:

1. source_data/example_inputs/ex_GVOGs.tsv

2. source_data/example_inputs/ex_GVEUKs.tsv

#### Run GV-host predictor:
```{bash}
python GVHoP.py --GVOGs_in=source_data/example_inputs/ex_GVOGs.tsv --GVEUKs_in=source_data/example_inputs/ex_GVEUKs.tsv --sample_ls=source_data/example_inputs/sample.ls --out_dir=temp_out/host_prediction/
```


## Results


## Contact

