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
wget "10.5281/zenodo.22199562" -O GVHoP_database_v1.0.zip
unzip GVHoP_database_v1.0.zip
```

## Example run
Activate GVHoP environment:

```{bash}
conda activate gvhop
```

Construct feature sets from sequence data:

```{bash}
python Extract/cmd_hmm.py \
    --hmm_dir /home/xinchang/hmm/bin/ \
    --in_dir=scratch/xinchang/girush02/girush02.26/source_data/run1/ \
    --db_dir=/scratch/xinchang/girush02/girush02.26/GVOGdb/GVHoP_GVOGs.hmm \
    --out_dir=/scratch/xinchang/girush02/girush02.26/hmmsearch/run1/ \
    --sh_dir=/scratch/xinchang/girush02/girush02.26/sh/hmmsearch/run1/ \
    --in_file_ext=faa \
    --out_file_ext=domout \
    --opt='hmmsearch --cpu 10 -E 1e-5 --domtblout' \
    --n_job=25
python3 /home/xinchang/pyscripts/execute_sh.py \
    --sh_dir=/scratch/xinchang/girush02/girush02.26/sh/hmmsearch/run1/ \
    --log_dir=/scratch/xinchang/girush02/girush02.26/execute/hmmsearch/run1/

# HGT signals
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

