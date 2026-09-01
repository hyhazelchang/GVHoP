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

## Example run
Activate GVHoP environment:

```{bash}
conda activate gvhop
```

1. Construct feature sets (Gene content) from sequence data:

1-1. This will generate the shell scripts for hmmsearch execution in sh_dir
```{bash}
python RunGVHoP/cmd_hmm.py --in_dir=../example_MAGs/ --db_dir=../GVHoP_database_v1.0/GVHoP_GVOGs.hmm --out_dir=../output_hmm/ --sh_dir=../sh/hmmsearch/ --in_file_ext=faa --out_file_ext=domout --opt='hmmsearch --cpu 10 -E 1e-5 --domtblout' --n_job=4
```

1-2. Execute hmmsearch
```{bash}
python RunGVHoP/execute_sh.py --sh_dir=../sh/hmmsearch/ --log_dir=../execute/hmmsearch/
```

1-3. Parse feature matrix
```{bash}
python RunGVHoP/parse_hmmsearch.py --in_dir=../output_hmm/ --in_file_ext=domout --out_dir=../source_data/example_inputs/ --out_preffix=ex_GVOGs --cogset=../source_data/features/GVHoP_GVOGs_all.tsv --seq_dir=../example_MAGs/ --seq_file_ext=faa
```

2. Construct feature sets (GV-euk signals) from sequence data:

2-1. This will generate the shell scripts for DIAMOND blastp execution in sh_dir
```{bash}
python3 RunGVHoP/cmd_diamond.py --task=blastp --in_dir=../example_MAGs/ --out_dir=../output_blastp/ --sh_dir=../sh/blastp/ --in_file_ext=faa --out_file_ext=txt --db_dir=../GVHoP_database_v1.0/GVHoP_GVEUKs.dmnd --outfmt=6 --threshold=8 --opt="--evalue 1e-5" --n_job=8
```

2-2. Execute diamond blastp
```{bash}
python RunGVHoP/execute_sh.py --sh_dir=../sh/blastp/ --log_dir=../execute/blastp/
```

2-3. Parse feature matrix
```{bash}
python RunGVHoP/parse_blasthits.py --blast_in=../output_blastp/ --in_file_ext=txt --column_names=../source_data/features/GVHoP_GVEUKs_all.tsv --outfile=../source_data/example_inputs/ex_GVEUKs.tsv
```


## Results


## Contact

