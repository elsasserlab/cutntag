import pandas as pd

from utils import (
    parse_stats_fields,
    parse_flagstat,
    parse_yaml,
    validate_config,
)


localrules:
    library_stats,
    stats_summary,
    spikein_factors,
    write_scaling_factor,
    unique_reads


TARGET = config.get("target")
SPIKEIN = config.get("spikein")
REF_LIBRARY = config.get("reflib")
QUALITY_CUTOFF = config.get("minquality", "20")
FRAG_LEN = 150

SAMPLES, = glob_wildcards("fastq/{sample}_R1.fastq.gz")

genomes_config = os.path.join(os.path.dirname(workflow.snakefile), "genomes.yaml")
genomes = parse_yaml(genomes_config)
validate_config(config, genomes)

multiqc_inputs = ([
        "reports/scaling_mqc.tsv",
        "reports/totals_stats_summary_mqc.tsv",
        "reports/pairs_stats_summary_mqc.tsv",
    ]
    + expand("stats/final/{sample}.insertsizes.txt", sample=SAMPLES)
    + expand("reports/fastqc/{sample}_R{read}_fastqc/fastqc_data.txt", sample=SAMPLES, read=(1, 2))
    + expand("log/bowtie/{sample}.target.all.bowtie.txt", sample=SAMPLES)
    + expand("log/bowtie/{sample}.spikein.all.bowtie.txt", sample=SAMPLES)
    + expand("stats/{sample}.dedup.target.metrics", sample=SAMPLES)
    + expand("stats/{sample}.dedup.spikein.metrics", sample=SAMPLES)
    + expand("final/peaks/{sample}_peaks.xls", sample=SAMPLES)
)

bigwigs = expand("final/bigwig/{sample}.bw", sample=SAMPLES) + \
          expand("final/bigwig/{sample}.unique.bw", sample=SAMPLES)

rpgc_bigwigs = expand("final/bigwig/{sample}.rpgc.bw", sample=SAMPLES) + \
          expand("final/bigwig/{sample}.unique.rpgc.bw", sample=SAMPLES)

bams_target = expand("final/target/{sample}.bam", sample=SAMPLES)
bams_spikein = expand("final/spikein/{sample}.bam", sample=SAMPLES)
bams = bams_target + bams_spikein

peaks = expand("final/peaks/{sample}_peaks.narrowPeak", sample=SAMPLES)


rule all:
    input:
        bams,
        bigwigs,
        peaks,
        "reports/multiqc_report.html"


rule no_spikein:
    input:
        bams_target,
        peaks,
        rpgc_bigwigs


rule no_bigwigs:
    input:
        bams,
        peaks,
        "reports/multiqc_report.html"


rule clean:
    shell:
        "rm -rf final/ tmp/ stats/ reports/ scaling/ log/"


rule multiqc:
    output: "reports/multiqc_report.html"
    input:
        multiqc_inputs,
        multiqc_config=os.path.join(os.path.dirname(workflow.snakefile), "multiqc.yaml")
    log:
        "log/multiqc_report.log"
    shell:
        "multiqc -f -o reports/ -c {input.multiqc_config} {input} 2> {log}"


rule fastqc_input:
    output:
        html="reports/fastqc/{sample}_fastqc.html",
        zip=temp("reports/fastqc/{sample}_fastqc.zip"),
        data="reports/fastqc/{sample}_fastqc/fastqc_data.txt",
    input:
        fastq="fastq/{sample}.fastq.gz"
    log:
        "log/1-fastqc/{sample}_fastqc.html.log"
    shell:
        "fastqc --extract -o reports/fastqc {input.fastq} > {log} 2>&1 "


# This removes unmapped reads directly, but keeps multimappers. We will remove
# multimappers from the counts for scaling, but keep them in place
# TODO: Maybe check for other bowtie settings or make them configurable.
rule align_target:
    threads: 8
    input:
        R1="fastq/{sample}_R1.fastq.gz",
        R2="fastq/{sample}_R2.fastq.gz",
    output:
        bam=temp("tmp/target/{sample}.all.bam")
    log:
        "log/bowtie/{sample}.target.all.bowtie.txt"
    params:
        reference=genomes[TARGET]["ref"],
    shell:
        "bowtie2"
        " -p {threads}"
        " -x {params.reference}"
        " -1 {input.R1}"
        " -2 {input.R2}"
        " "
        " --fast"
        " 2> {log}"
        " "
        "| samtools sort -@ {threads} "
        "| samtools view -F 4 -b - > {output.bam}"


rule align_spikein:
    threads: 8
    input:
        R1="fastq/{sample}_R1.fastq.gz",
        R2="fastq/{sample}_R2.fastq.gz",
    output:
        bam=temp("tmp/spikein/{sample}.all.bam")
    log:
        "log/bowtie/{sample}.spikein.all.bowtie.txt"
    params:
        reference=genomes[SPIKEIN]["ref"],
    shell:
        "bowtie2"
        " -p {threads}"
        " -x {params.reference}"
        " -1 {input.R1}"
        " -2 {input.R2}"
        " --fast"
        " 2> {log}"
        " "
        "| samtools sort -@ {threads} "
        "| samtools view -F 4 -b - > {output.bam}"


rule mark_duplicates_target:
    input:
        bam="tmp/target/{sample}.all.bam"
    output:
        dedup_bam=temp("tmp/target/{sample}.dedup.bam"),
        dedup_metrics="stats/{sample}.dedup.target.metrics"
    log:
        "log/target/{sample}.markdups.txt.log"
    shell:
        "picard MarkDuplicates"
        " REMOVE_DUPLICATES=TRUE"
        " I={input.bam}"
        " O={output.dedup_bam}"
        " M={output.dedup_metrics}"
        " 2> {log}"


rule mark_duplicates_spikein:
    input:
        bam="tmp/spikein/{sample}.all.bam"
    output:
        dedup_bam=temp("tmp/spikein/{sample}.dedup.bam"),
        dedup_metrics="stats/{sample}.dedup.spikein.metrics"
    log:
        "log/target/{sample}.markdups.txt.log"
    shell:
        "picard MarkDuplicates"
        " REMOVE_DUPLICATES=TRUE"
        " I={input.bam}"
        " O={output.dedup_bam}"
        " M={output.dedup_metrics}"
        " 2> {log}"


rule remove_exclude_regions_target:
    threads: 8
    input:
        bam="tmp/target/{sample}.dedup.bam",
    params:
        exclude=genomes[TARGET]["exclude"],
    output:
        bam="final/target/{sample}.bam"
    shell:
        "bedtools"
        " intersect"
        " -v"
        " -abam {input.bam}"
        " -b {params.exclude}"
        " > {output.bam}"


rule remove_exclude_regions_spikein:
    threads: 8
    input:
        bam="tmp/spikein/{sample}.dedup.bam",
    params:
        exclude=genomes[SPIKEIN]["exclude"],
    output:
        bam="final/spikein/{sample}.bam"
    shell:
        "bedtools"
        " intersect"
        " -v"
        " -abam {input.bam}"
        " -b {params.exclude}"
        " > {output.bam}"


rule unique_reads:
    input:
        bam="final/{sample}.bam"
    output:
        bam="final/{sample}.unique.bam"
    params:
        quality=QUALITY_CUTOFF
    shell:
        "samtools view -bh -q {params.quality} {input.bam} > {output.bam}"


rule samtools_index:
    threads: 4
    input:
        bam="{sample}.bam"
    output:
        index=temp("{sample}.bai")
    shell:
        "samtools index -@ {threads} -b {input.bam} {output.index}"


rule samtools_flagstat:
    threads: 4
    input:
        bam="{sample}.bam",
        index="{sample}.bai"
    output:
        stats=temp("{sample}.flagstat.txt")
    shell:
        "samtools flagstat -@ {threads} {input.bam} > {output.stats}"


rule library_stats:
    output:
        total=temp("stats/{sample}.total.txt"),
        pairs=temp("stats/{sample}.pairs.txt"),
    input:
        target_mapped="tmp/target/{sample}.all.flagstat.txt",
        target_dedup="tmp/target/{sample}.dedup.flagstat.txt",
        target_fltd="final/target/{sample}.flagstat.txt",
        target_fltd_unique="final/target/{sample}.unique.flagstat.txt",
        spikein_mapped="tmp/spikein/{sample}.all.flagstat.txt",
        spikein_dedup="tmp/spikein/{sample}.dedup.flagstat.txt",
        spikein_fltd="final/spikein/{sample}.flagstat.txt",
        spikein_fltd_unique="final/spikein/{sample}.unique.flagstat.txt",
    run:
        d_total = dict(library=f"{wildcards.sample}")
        d_pairs = dict(library=f"{wildcards.sample}")

        for flagstat, name in [
            (input.target_mapped, "target_mapped"),
            (input.target_dedup, "target_dedup"),
            (input.target_fltd, "target_fltd"),
            (input.target_fltd_unique, "target_fltd_unique"),
            (input.spikein_mapped, "spikein_mapped"),
            (input.spikein_dedup, "spikein_dedup"),
            (input.spikein_fltd, "spikein_fltd"),
            (input.spikein_fltd_unique, "spikein_fltd_unique")
        ]:
            d_total[name] = parse_flagstat(flagstat)["mapped"]
            d_pairs[name] = parse_flagstat(flagstat)["pairs"]

        with open(output.total, "w") as f:
            print(*d_total.keys(), sep="\t", file=f)
            print(*d_total.values(), sep="\t", file=f)

        with open(output.pairs, "w") as f:
            print(*d_pairs.keys(), sep="\t", file=f)
            print(*d_pairs.values(), sep="\t", file=f)


rule insert_size_metrics:
    threads: 4
    output:
        txt="stats/final/{name}.insertsizes.txt",
        pdf="stats/final/{name}.insertsizes.pdf",
    input:
        bam="final/target/{name}.bam"
    log:
        "log/final/{name}.insertsizes.txt.log"
    shell:
        "picard"
        " CollectInsertSizeMetrics"
        " I={input.bam}"
        " O={output.txt}"
        " HISTOGRAM_FILE={output.pdf}"
        " MINIMUM_PCT=0.5"
        " STOP_AFTER=10000000"
        " 2> {log}"


rule stats_summary:
    output:
        totals="reports/totals_stats_summary_mqc.tsv",
        pairs="reports/pairs_stats_summary_mqc.tsv"
    input:
        totals=expand("stats/{sample}.total.txt", sample=SAMPLES),
        pairs=expand("stats/{sample}.pairs.txt", sample=SAMPLES)
    run:
        header = [
            "library",
            "target_mapped",
            "target_dedup",
            "target_fltd",
            "target_fltd_unique",
            "spikein_mapped",
            "spikein_dedup",
            "spikein_fltd",
            "spikein_fltd_unique",
        ]

        with open(output.totals, "w") as f:
            print(*header, sep="\t", file=f)
            for stats_file in sorted(input.totals):
                summary = parse_stats_fields(stats_file)
                row = [summary[k] for k in header]
                print(*row, sep="\t", file=f)

        with open(output.pairs, "w") as f:
            print(*header, sep="\t", file=f)
            for stats_file in sorted(input.pairs):
                summary = parse_stats_fields(stats_file)
                row = [summary[k] for k in header]
                print(*row, sep="\t", file=f)


rule spikein_factors:
    input:
        totals="reports/pairs_stats_summary_mqc.tsv",
    output:
        report="reports/scaling_mqc.tsv",
    params:
        genome_size=genomes[TARGET]["size"],
        fraglen=FRAG_LEN,
        reflib=REF_LIBRARY
    run:
        use_column = "fltd_unique"
        df = pd.read_table(input.totals)
        df.set_index("library", inplace = True, drop = False)
        min_spikein = min(df[f"spikein_{use_column}"])
        df["spikein_inv"] = min_spikein / df[f"spikein_{use_column}"]
        df["spikein_raw_ratio"] = df[f"target_{use_column}"] / df[f"spikein_{use_column}"]
        df["scaling_factor"] = (float(params.genome_size) * df["spikein_raw_ratio"]) / (params.fraglen*df[f"target_{use_column}"]) # changed spikein_factor to spike_in_raw_ratio, SE 2023-07-09
        df["scaling_factor_rpgc"] = df["spikein_raw_ratio"] / df.loc[params.reflib]["spikein_raw_ratio"]
        df.to_csv(output.report, sep="\t", index = False)


rule write_scaling_factor:
    input:
        report="reports/scaling_mqc.tsv"
    output:
        factor="scaling/{sample}.factor",
        factor_unique="scaling/{sample}.unique.factor"
    run:
        df = pd.read_table(input.report)
        df.set_index("library", inplace = True)
        n = df.loc[wildcards.sample]["scaling_factor_rpgc"]
        with open(output.factor, "w") as f:
            f.write(f"{n}")
        with open(output.factor_unique, "w") as f:
            f.write(f"{n}")


rule generate_scaled_bigwig:
    threads: 8
    input:
        factor="scaling/{sample}.factor",
        bam="final/target/{sample}.bam",
        index="final/target/{sample}.bai"
    output:
        bigwig="final/bigwig/{sample}.bw"
    log:
        "log/final/{sample}.bw.log"
    params:
        genome_size=genomes[TARGET]["size"],
        fraglen=FRAG_LEN
    shell:
        "bamCoverage"
        " -p {threads}"
        " -b {input.bam}"
        " --binSize 1"
        " --normalizeUsing RPGC"
        " --effectiveGenomeSize {params.genome_size}"
        " --extendReads {params.fraglen}"
        " --scaleFactor $(< {input.factor})"
        " -o {output.bigwig}"
        " 2> {log}"


rule generate_unscaled_bigwig:
    threads: 8
    input:
        bam="final/target/{sample}.bam",
        index="final/target/{sample}.bai"
    output:
        bigwig="final/bigwig/{sample}.rpgc.bw"
    log:
        "log/final/{sample}.bw.log"
    params:
        genome_size=genomes[TARGET]["size"],
        fraglen=FRAG_LEN
    shell:
        "bamCoverage"
        " -p {threads}"
        " -b {input.bam}"
        " --binSize 1"
        " --normalizeUsing RPGC"
        " --effectiveGenomeSize {params.genome_size}"
        " --extendReads {params.fraglen}"
        " -o {output.bigwig}"
        " 2> {log}"


# TODO: Add option for --broad --broad-cutoff 0.05 or whatever value.
rule call_peaks:
    threads: 4
    input:
        bam="final/target/{sample}.bam",
        bai="final/target/{sample}.bai"
    output:
        peaks="final/peaks/{sample}_peaks.narrowPeak",
        peaks_xls="final/peaks/{sample}_peaks.xls"
    params:
        genome_size=genomes[TARGET]["size"],
        peaksdir="final/peaks/"
    shell:
        "macs3 callpeak"
        " -t {input.bam}"
        " -n {wildcards.sample}"
        " -s {params.genome_size}"
        " -f BAMPE"
        " --outdir {params.peaksdir}"
