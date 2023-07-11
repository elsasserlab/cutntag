import os
import yaml

def parse_stats_fields(stats_file):
    """
    Parse contents of a stats file created in the stats rule, which consists
    of a header line and a values line.

    Return a dictionary where keys are defined by the values of the header.
    """
    with open(stats_file) as f:
        header = f.readline().strip().split("\t")
        values = f.readline().strip().split("\t")
    result = {key.lower(): value for key, value in zip(header, values)}
    return result


def parse_flagstat(path):
    counts = {}
    with open(path) as f:
        for line in f:
            if " in total" in line:
                counts["total"] = int(line.split(maxsplit=1)[0])
            elif " mapped (" in line:
                counts["mapped"] = int(line.split(maxsplit=1)[0])
            elif " primary duplicates" in line:
                counts["duplicates"] = int(line.split(maxsplit=1)[0])
            elif " properly paired" in line:
                counts["pairs"] = int(line.split(maxsplit=1)[0])

    return counts


def parse_yaml(path):
    with open(path, "r") as fi:
        return yaml.safe_load(fi)


def check_parameter_exists(config, key):
    if key not in config.keys():
        msg = f"Missing required parameter {key}."
        raise ValueError(msg)


def validate_genome_size(value):
    try:
        gsize = int(value)
        if gsize <= 0:
            msg = "Genome size must be a positive value"
            raise ValueError(msg)

    except ValueError:
        msg = "An integer value is required as effective genome size"
        raise ValueError(msg)


def validate_reference(ref):
    refbase = ref["ref"]
    bowtie_file = f"{refbase}.1.bt2"
    if not os.path.isfile(bowtie_file):
        msg = f"Bowtie2 index file {bowtie_file} not found for reference."
        raise FileNotFoundError(msg)

    bedfile = ref["exclude"]
    if not os.path.isfile(bedfile):
        msg = f"Exclude BED file {bedfile} not found for reference."
        raise FileNotFoundError(msg)


    validate_genome_size(ref["size"])


def validate_config(config, genomes):
    check_parameter_exists(config, "target")
    
    target = config.get("target")
    if target not in genomes.keys():
        msg = f"Unknown target reference {target}. Please add it to the genomes list"
        raise ValueError(msg)
    validate_reference(genomes[target])

    spikein = config.get("spikein")
    if spikein:
        if spikein not in genomes.keys():
            msg = f"Unknown spikein reference {spikein}. Please add it to the genomes list"
            raise ValueError(msg)
        validate_reference(genomes[spikein])

    