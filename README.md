# Single molecule methylation profiles of cell-free DNA in cancer with nanopore sequencing
Code repository of scripts used in Lau, et al. 2022. Specifically, scripts used to generate in silico admixtures, as well as scripts for single molecule methylation classification of reads against reference profiles are included here. Some R scripts require data munging of BAM files (e.g. extraction of read lengths); these are not included but the expected formatting is described.

## Contents
_scripts use standard functions and are likely to work with other package versions if using Python 3+ and R 4+)_

**nanopore_mix_bams.py** - python script for mixing bam files together using a random seed and proportion. The ground truth source of each read is also record in a log file. This has been tested with the following:
  - python 3.8.12 (base packages argparse, csv, sys, random, os)
  - pysam 0.18.0
  ```
  usage: nanopore_mix_bams.py [-h] --bam1 BAM1 --bam2 BAM2 --out_prefix OUT_PREFIX [--frac FRAC] [--nreads NREADS] [--threads THREADS] [--seed SEED]

  optional arguments:
  -h, --help            show this help message and exit
  --bam1 BAM1           first bam file
  --bam2 BAM2           second bam file
  --out_prefix OUT_PREFIX
                        output prefix
  --frac FRAC           mixing fraction (p <= bam 1/total)
  --nreads NREADS       number of reads to generate
  --threads THREADS     number of threads
  --seed SEED           random seed
  ```

**nanopore_dump_reads.py** - python script to dump a bam file with modified base tags into a flat file. This has been tested with the following:
  - python 3.8.12 (base packages argparse, csv, sys)
  - pandas 1.3.4
  - ModBam 0.3.3 and modbampy 0.3.2 ([ModBam/modbampy is maintained by Oxford Nanopore](https://github.com/epi2me-labs/modbam2bed))

  ```
  usage: nanopore_dump_reads.py [-h] --bam BAM --fasta FASTA --out OUT [--nthreads NTHREADS]

  optional arguments:
  -h, --help           show this help message and exit
  --bam BAM            bam file
  --fasta FASTA        fasta file
  --out OUT            output table
  --nthreads NTHREADS  number of threads
  ```

**nanopore_classify_reads.R** - R notebook-style script used in Rstudio to classify reads from nanopore_dump_reads.py against reference methylomes. This has been tested with the following:
  - Rstudio Server 1.4.1717
  - R 4.1.0
  - data.table 1.14
  - dplyr 1.0.7
  - ggplot2 3.3.5
  - readxl 1.3.1
  - tidyverse 1.3.1

**nanopore_gene_analysis.R** - R notebook-style script used in Rstudio to perform gene-level methylation anaysis and plotting. This has been tested with the following:
  - Rstudio Server 1.4.1717
  - R 4.1.0
  - data.table 1.14
  - dplyr 1.0.7
  - ggplot2 3.3.5
  - readxl 1.3.1
  - tidyverse 1.3.1
  - pheatmap 1.0.12

**nanopore_size_analysis.R** - R notebook-style script used in Rstudio to compute insert size distributions of cfDNA nanopore sequencing data. This has been tested with the following:
  - Rstudio Server 1.4.1717
  - R 4.1.0
  - data.table 1.14
  - dplyr 1.0.7
  - ggplot2 3.3.5
  - tidyverse 1.3.1
