# TRIM Identification

---

## Authors

Sophie Maiwald  
Ferdinand Maiwald  
Tony Heitkam

---

## TRIM Identification Introduction

This repository provides a detailed workflow for the detection and annotation of **non-autonomous LTR retrotransposons (TRIMs)** 

The overall workflow consists of three major steps illustrated in the figure below.

![TRIM identification workflow](https://github.com/TRIMpossible/TRIM_identification/blob/main/images/TRIM_identification_workflow.png)

---

## Requirements

### Software

The following software is required:

- Python ≥ 3.8
- BLAST+
- LTR_Finder
- Flexidot
- MCL or SiLix (for clustering)

### Python dependencies

Install required packages:

```bash
pip install biopython bcbio-gff
```

---

## Step 1: *De novo* TRIM identification 

We recommend running LTR_Finder within a conda environment.

For installation instructions please refer to Bioconda: [LTR_Finder instructions](https://bioconda.github.io/recipes/ltr_finder/README.html)

LTR_Finder parameters must be adjusted to optimize TRIM identification:
- reduce minimal LTR length (-l)
- reduce minimal element distance (-d)
- disable detection of coding domains (-F)

For PBS prediction we recommend using plant tRNAs from: [tRNA database](http://seve.ibmp.unistra.fr/plantrna)

```bash
conda activate ltr_finder

ltr_finder \
 -w 2 \
 -d 30 \  -D 10000 \
 -l 30 \  -L 5000 \
 -r 10 \ -s tRNA_database.fasta \
 -E \  -F 00001110000 \
 genome_assembly.fasta \  > LTR_Finder_out.txt
```

LTR_Finder output can be extracted as nucleotide sequences using:

```bash
helper_scripts/extract_LTR_Finder_fasta.py  \
  --ltr LTR_Finder_out.txt \ --genome genome_assembly.fasta \
  --output LTR_Finder_output.fasta
```

### LTR_Finder output filtering 

To estimate full copy numbers of non-autonomous elements, LTR_Finder outputs must be further processed.

*Remove sequences with coding domains*

Perform a BLASTx search of LTR_Finder outputs against coding domains of autonomous LTR retrotransposons.

BLAST+ installation instructions: [BLAST install](https://www.ncbi.nlm.nih.gov/books/NBK569861/)

database of LTR-retrotransposon domains: [REXdb](http://repeatexplorer.org/?page_id=918)

```bash
makeblastdb \
 -in protein_domain_database.fasta \ -parse_seqids \
 -dbtype prot \ -out database_name

blastx \
 -query LTR_Finder_output.fasta \
 -db database_name \
 -outfmt 6 \ -out results.m6
```

*Exclude coding-domain hits*

Exclude all LTR_Finder output sequences with similarity to coding domains.

*Flexidot analysis*

Run Flexidot on BLAST-filtered LTR_Finder output.

Installation instructions: [flexidot install](https://github.com/flexidot-bio/flexidot)

```bash
conda activate flexidot

flexidot -f 1 -k 7 -p 0 -r N -E 7
```
We recommend using the Flexidot classifier: [flexidot classifier](https://github.com/tudipffmgt/Flexidot-Classifier)

This classifier estimates full-length TRIM sequences and reduces manual curation effort.

Filtered TRIM sequences can then be used for genome-wide quantification (Step 2).

---

## Step 2: Quantification and full length annotation

Genome-wide sequence searches can be performed using BLAST+.

```bash
makeblastdb \
 -in TRIM_reference_for_quantification.fasta \ -parse_seqids \
 -dbtype nucl \
 -out database_quantification

blastn \
 -task megablast \
 -query genome_assembly.fasta \ -db database_quantification \
 -word_size 12 \
 -outfmt 6 \ -out results_quantification.m6
```

### Quantification output processing 

For quantification output processing we recommend using scripts stored under: ```bash blast_output_processing/ ```

These scripts perform:
- Filtering BLAST output by query
- Conversion of BLAST output to GFF
- Merging of nearby hits (default 5 kb distance)
- Extraction of FASTA sequences from GFF coordinates
- Optional extraction with flanking regions
- Structural splitting of sequences
- Self-BLAST validation of candidate TRIM elements

### Run workflow

The complete workflow can be executed using:

```bash
python quantification_output_processing.py \
 --blast results_quantification.m6 \ 
 --genome genome.fasta \
 --workdir run1 \ 
 --run-selfblast \ --selfblast-input split \ --selfblast-combined-gff
```

This will execute: BLAST output → GFF → merged GFF → FASTA → self-BLAST

### Detailed Annotation

More detailed and accurate full-length annotation can be performed using Flexidot and the Flexidot classifier again.

---

## Step 3: TRIM sequence clustering

Clustering of annotated full-length sequences is best achieved using an all-against-all BLAST search.

```bash
makeblastdb \
 -in TRIM_full_lengths.fasta \
 -parse_seqids \
 -dbtype nucl \ -out TRIM_allVall

blastn \
 -query TRIM_full_lengths.fasta \
 -db TRIM_allVall \
 -outfmt 6 \ -out TRIM_allVall_blast.m6
```

Clustering can be performed using either:

- *SiLix* [Silix](https://lbbe.univ-lyon1.fr/fr/SiLix)

-  *Markov Clustering (MCL)* [MCL](https://anaconda.org/bioconda/mcl)

```bash
conda install -c bioconda mcl
```

---
 
## Further helpfull scripts

*Filtering duplette sequences after clustering*

This script will elongate full length sequences with flanking regions and performs MAFFT alignments in order to filter duplicated hits: 

```bash
python scanning_cluster_for_duplicates.py
```

*Annotation of target site duplications (TSD)*

Annotation of TSDs will be performed on full length sequences with flanking regions.

In order for the script to work a gff file is needed providing an "element" annotation from 5'LTR start to 3' LTR end. 

```bash
python find_tsds.py --gff_input --fasta_input --output_prefix)
```
