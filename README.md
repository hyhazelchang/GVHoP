# GVHoP: Giant Virus Host Predictor

GVHoP predicts eukaryotic host groups for giant virus metagenome-assembled genomes (GVMAGs) from their predicted proteins. It combines protein-family HMM matches and similarities to eukaryotic proteins using trained XGBoost classifiers and a neural meta-model.

## Installation

Clone the repository and create the environment described in `env.yml`:

```bash
git clone https://github.com/hyhazelchang/GVHoP.git
cd GVHoP
conda env create -f env.yml
conda activate gvhop
chmod u+x gvhop
```

The environment includes Python 3.11, HMMER, DIAMOND, XGBoost, scikit-learn, and CPU-only PyTorch. A GPU is not required. The repository includes the trained models in `XGBclf/` and `NNclf/`; keep these directories and `source_data/` beside the `gvhop` executable. The checkout occupies about 7.1 GiB, including Git history.

Download the database archive from the [official Zenodo record](https://zenodo.org/records/22199562), verify its published checksum, and extract it:

```bash
curl --fail --location \
  https://zenodo.org/api/records/22199562/files/GVHoP_database_v1.0.zip/content \
  --output GVHoP_database_v1.0.zip
printf '%s\n' 'f5ab380ffa252de59d1d00f671a887f6  GVHoP_database_v1.0.zip' | md5sum --check
unzip GVHoP_database_v1.0.zip
```

The archive is 205,381,301 bytes and extracts the two required files into `db/`: `GVHoP_GVOGs.hmm` and `GVHoP_GVEUKs.dmnd`. Their combined size is 2,431,490,723 bytes. Allow space for the software environment and run outputs in addition to the checkout and database.

## Inputs and execution

Supply one protein FASTA file per genome. GVHoP does not predict genes or translate nucleotide FASTA files. List the genomes to process in a plain-text sample file, with one sample identifier per line and no header or file extension. For example, sample `virus_A` and extension `faa` identify `virus_A.faa` in the input directory.

Sample identifiers must be unique and nonempty. They cannot contain `/`, `\`, or NUL, or equal `.` or `..`. Interior spaces are allowed. Only files named in the sample list are searched; each must exist and be nonempty.

Run the supplied protein examples from the repository directory:

```bash
./gvhop \
  --in_dir example_MAGs \
  --in_file_ext faa \
  --db_dir db \
  --sample_ls source_data/example_inputs/sample.ls \
  --out_dir gvhop_out \
  --cpu 8
```

| Argument | Meaning |
| --- | --- |
| `--in_dir` | Protein FASTA directory; defaults to `example_MAGs` beside the executable. |
| `--in_file_ext` | Required filename extension, such as `faa`, without a leading dot. |
| `--db_dir` | Required directory containing both nonempty database files. |
| `--sample_ls` | Required sample-list file. |
| `--out_dir` | Required output directory; it must be absent or empty. |
| `--cpu` | Required positive integer specifying the total CPU budget. |

Paths can be relative or absolute and do not need trailing slashes. Quote paths containing spaces. You can invoke `gvhop` from another working directory: its model and feature files are resolved relative to the executable, while relative paths passed as arguments are resolved from your working directory.

Searches run sequentially, with one external search process at a time. `--cpu` bounds the search workers and the XGBoost, PyTorch, and native numerical-library threads. HMMER reserves one thread from that budget for coordination. The supplied environment runs inference on CPUs; GPU use requires compatible libraries and hardware.

Use a fresh output directory for each run. GVHoP rejects nonempty output directories to prevent mixing results from different inputs or failed attempts. A failed external search stops the pipeline with a nonzero exit status. Inspect its per-sample log before retrying in a new directory. A completed search with no hits is valid: its features are zero, and inference still runs.

## Outputs

Feature and prediction tables use tabs and identify samples in the `testset` column.

| Path inside the output directory | Contents |
| --- | --- |
| `FeatureSets/GVOGs_all.tsv` | Counts for 8,293 protein-family features. |
| `FeatureSets/GVEUKs_all.tsv` | Similarity scores for 57,250 eukaryotic protein features. |
| `FeatureSets/GVOGs_top.tsv` | The 444 protein-family features also passed to the meta-model. |
| `FeatureSets/GVEUKs_top.tsv` | The 447 similarity features also passed to the meta-model. |
| `Predictions/prob_all.tsv` | `testset` followed by probabilities for 21 host labels. |
| `Predictions/h_prob_all.tsv` | `testset`, `pred_label`, `pred_hierarchy`, `level0_prob`, `level1_prob`, `level2_prob`. |
| `Predictions/pred_out_level0.tsv` | One row per sample; accepted prediction counts for the 21 finest host labels. |
| `Predictions/pred_out_level1.tsv` | Accepted prediction counts for 11 intermediate host groups. |
| `Predictions/pred_out_level2.tsv` | Accepted prediction counts for the five broadest host groups. |
| `hmmsearch/<sample>.domout` | HMMER domain-search results. |
| `hmmsearch/<sample>.domout.parsed` | Best-hit records with `protein_id`, `accession`, `best_hit`, `aln_start`, `aln_end`, `aln_length`, `score`, and `type`. |
| `diamond_blastp/<sample>.txt` | Headerless DIAMOND results in its 12-column format 6. |
| `execute/hmmer/<sample>.log` | HMMER command output for one sample. |
| `execute/diamond_blastp/<sample>.log` | DIAMOND command output for one sample. |

Each of the 100 model iterations combines six XGBoost classifiers with the same neural meta-model. Each probability table has 100 rows per sample, without an iteration identifier. The 21 probability columns in `prob_all.tsv` are ordered as follows; this file has no `pred_label` column:

```text
Aves, Bathycoccaceae, Blastocladiomycetes, Chlorellaceae, Choanocafe,
Chytridiomycetes, Coccolithophyceae, Discoba, Discosea, Insecta,
Malacostraca, Mamiellaceae, Mammalia, Otherchlorophyta, Otherhaptophyta,
Otherinvertebrate, Othervertebrate, Phaeophyceae, Prymnesiaceae,
Stramenopile, Tubulinea
```

These probabilities sum to one per row, subject to floating-point rounding. They describe the model's 21 host labels, which use mixed taxonomic ranks.

In `h_prob_all.tsv`, `pred_hierarchy` contains three comma-separated host names, ordered from finest (`level0`) to broadest (`level2`). `pred_label` contains three corresponding zero-based indices into the alphabetically sorted labels at each level. Starting with the highest-probability finest label, GVHoP follows its branch toward broader groups. It accepts a level when the summed probability of that group's finest labels is at least 0.75. A rejected level uses `-` in both hierarchy fields and its probability field. The broadest groups are `Algae`, `Amoeba`, `Fungi`, `Heteroflagellate`, and `Metazoa`.

Illustrative example: `-,Amoebozoa,Amoeba` with probabilities `-`, `0.80`, and `0.80` accepts the two broader groups but abstains at the finest level. These values are not a measured result for a supplied genome.

The three count tables aggregate accepted labels across the 100 predictions. Their taxon columns contain integer counts from 0 to 100. A row can sum to less than 100 because abstentions are omitted; an all-zero row means no prediction passed the threshold at that level. Missing votes and zero-hit searches do not by themselves indicate a failed run.

## Regression tests

Run the parser, input-validation, path-handling, failure-propagation, and vote-counting tests in the installed GVHoP environment:

```bash
python -m unittest discover -s tests -v
```

The tests use the standard-library `unittest` runner and the application dependencies. They do not require the trained models or downloaded databases.

## Contact and citation

For questions, contact [Hsin-Ying Chang](https://github.com/hyhazelchang).

Machine learning prediction of eukaryotic hosts for giant viruses. [https://doi.org/10.64898/2026.09.08.750034](https://doi.org/10.64898/2026.09.08.750034)
