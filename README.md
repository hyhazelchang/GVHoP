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

Construct feature sets (Gene content) from sequence data:
This will generate the shell scripts for hmmsearch execution in sh_dir

```{bash}
python RunGVHoP/cmd_hmm.py --in_dir=../example_MAGs/ --db_dir=../GVHoP_database_v1.0/GVHoP_GVOGs.hmm --out_dir=../output_hmm/ --sh_dir=../sh/hmmsearch/ --in_file_ext=faa --out_file_ext=domout --opt='hmmsearch --cpu 10 -E 1e-5 --domtblout' --n_job=4
```

Execute hmmsearch
```{bash}
python RunGVHoP/execute_sh.py --sh_dir=../sh/hmmsearch/ --log_dir=../execute/hmmsearch/
```

## Parse feature matrix
python RunGVHoP/parse_hmmsearch.py --in_dir=../output_hmm/ --in_file_ext=domout --out_dir=../source_data/example_inputs/ --out_preffix=ex_GVOGs --cogset=../source_data/features/GVHoP_GVOGs_all.tsv --seq_dir=../example_MAGs/ --seq_file_ext=faa


# GV-euk signals
python3 /home/xinchang/pyscripts/cmd_diamond.py \
    --diamond_dir=/home/xinchang/softwares/diamond \
    --task=blastp \
    --in_dir=/scratch/xinchang/girush02/girush02.26/source_data/run1/ \
    --out_dir=/scratch/xinchang/girush02/girush02.26/blast_eukHGT_v3/run1/ \
    --sh_dir=/scratch/xinchang/girush02/girush02.26/sh/blast_eukHGT_v3/run1/ \
    --in_file_ext=faa \
    --out_file_ext=txt \
    --db_dir=/scratch/xinchang/girush02/girush02.26/GVEUKdb/GVHoP_GVEUKs \
    --outfmt=6 \
    --threshold=10 \
    --opt="--evalue 1e-5" \
    --n_job=10;
nice -n 5 python3 /home/xinchang/pyscripts/execute_sh.py \
    --sh_dir=/scratch/xinchang/girush02/girush02.26/sh/blast_eukHGT_v3/run1/ \
    --log_dir=/scratch/xinchang/girush02/girush02.26/execute/blast_eukHGT_v3/run1/

```



## Results

## Contact

