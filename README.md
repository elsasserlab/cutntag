# Simple CUT&Tag workflow

A simple Snakemake workflow to process paired-end CUT&amp;Tag libraries and produce spike-in scaled bigWig files.

## Workflow description

This is a very simple workflow that performs the following steps:

1.  Align with `bowtie2` to reference and spikein reference genomes.
2.  Deduplicate with `Picard`.
3.  Filter regions from an exclude loci set (BED file) using `bedtools intersect`.
4.  Collect insert size metrics of these alignments with `Picard`.
5.  Filter for uniquely mapping reads (MQ over certain threshold - default \> 20) using `samtools view`.
6.  Count alignments both total and proper pairs, making two output tables using `samtools flagstat`.
7.  Calculate spikein factors (see scaling below).
8.  Make scaled bigWig files using spikein factors using `deepTools bamCoverage`.
9.  Call peaks with `MACS3` on each sample (MACS3 has a default non-background method now for ATAC seq that should
    work with CUT&Tag better than the ChIP-seq specific one).

## Spikein scaling

This workflow produces bigWig files where the sample reference is scaled to RPGC 1x genome coverage using deepTools.
In order to achieve this, the spike-in factor is corrected to account for this. deepTools is called
with `--effectiveGenomeSize` set accordingly to the target genome, and a `--scaleFactor` parameter that is calculated
as the ratio: `target_reads / spikein_reads`, normalized to the reference ratio
`ref_target_reads / ref_spikein_reads` (proper pair counts, after deduplication and exclude regions, to
avoid competition for possibly shared repetitive regions).

Note that for the reference this second factor always equals to 1, hence, reference library has always global
mean of 1. The rest of the libraries are relative to it.

## Setup

1. Create a conda environment using the `environment.yml` file provided.

```
mamba env create -n cutntag -f environment.yml
```

2. Activate the environment

```
conda activate cutntag
``` 

3. Copy the `cutntag.yaml` file to your running experiment directory (where the `fastq` directory with the FASTQ files is contained).
4. Run snakemake defining target, spikein and reflib parameters:

```
snakemake -p -s /path/to/cutntag/Snakemake --rerun-incomplete --config target=mm39 spikein=dm6 reflib=df_library
```

Note that `target`, `spikein` and `reflib` parameters can also be defined in the `cutntag.yaml` file, and that their values
must match one of the references in the same file. For an example of those files, check the 

### Configuration

The workflow will run on all the FASTQ file pairs provided in the `fastq` directory in the running directory.

Additionally, it needs a cutntag.yaml file that contains paths to the reference genome files, bowtie indexes and exclude BED files.
It takes as parameters `target` and `spikein` which should match a genome in the `config.yaml` file, and `reflib` as the
prefix of the fastq file pair that serves as spikein reference. For instance, for `tf_library_R1.fastq.gz`, `tf_library_R2.fastq.gz`,
the value would be `tf_library` (without the `_R1.fastq.gz`).




