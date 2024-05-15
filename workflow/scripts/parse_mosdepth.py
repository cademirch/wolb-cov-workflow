import pandas as pd
from pathlib import Path


def parse_flagstat(flagstat_file: Path) -> dict:
    mapped_reads = None
    duplicates = None
    total_reads = None

    with open(flagstat_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            if "total" in line:

                total_reads = int(line.split()[0])
            elif "mapped" in line:

                mapped_reads = int(line.split()[0])
            elif "duplicates" in line:

                duplicates = int(line.split()[0])

            # Break the loop if all values are found
            if all([x is not None for x in (mapped_reads, duplicates, total_reads)]):
                break

    return {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "duplicate_reads": duplicates,
    }


def parse_mosdepth(file, groups, roles):
    """
    groups = {"dmel":[<dmel autosome>], ...}
    groups should only be the genomes for the given sample
    """
    inverse_groups = {}
    out = {}
    for group, contigs in groups.items():
        out[group] = 0
        for contig in contigs:
            inverse_groups[contig] = group
    with file.open() as f:
        for line in f:
            line = line.strip().split()
            chrom = line[0]
            if len(chrom.split("_region")) == 2:
                chrom = chrom.split("_region")[0]
                cov = float(line[3])
                if chrom in inverse_groups:
                    group = inverse_groups[chrom]
                    out[group] += cov
    # Divide each group by number of contigs
    for group, contigs in groups.items():
        out[f"{group}_mean_depth"] = out.pop(group) / len(contigs)  # pop to rename key

    if len(groups) > 1:
        host = [g for g in groups if g in roles["hosts"]][0]
        symbs = [g for g in groups if g in roles["symbionts"]]
        for symb in symbs:
            try:
                out[f"{symb}_titer"] = (
                    out[f"{symb}_mean_depth"] / out[f"{host}_mean_depth"]
                )
            except ZeroDivisionError:
                out[f"{symb}_titer"] = 0
    return out


def main():
    groups = snakemake.params["groups"]  # noqa: F821
    roles = snakemake.params["roles"]  # noqa: F821
    depths = snakemake.input.depths  # noqa: F821
    flagstat_files = snakemake.input.dedup_stats  # noqa: F821
    sample_sheet = snakemake.params.sample_sheet  # noqa: F821
    outfile = Path(snakemake.output["csv"])  # noqa: F821
    samples = pd.read_json(sample_sheet, orient="records")
    sampleids = samples["SampleID"].unique().tolist()

    file_dict = {}

    for sampid in sampleids:
        depth = [f for f in depths if sampid in f][0]
        flag = [f for f in flagstat_files if sampid in f][0]
        file_dict[sampid] = (depth, flag)

    dicts = []

    for sample_id, files in file_dict.items():

        depth_file, flagstat_file = files
        row = samples[samples["SampleID"] == sample_id].iloc[0]
        host = row["host"]
        infection = row["infection"]
        genomes = []
        if infection:
            genomes.extend(infection)

        genomes.append(host)
        print(sample_id, host, infection, genomes)
        genome_groups = dict((k, groups[k]) for k in genomes)

        out = {"SampleID": sample_id, "sequencing_run": row["sequencing_order_id"]}
        out.update(parse_mosdepth(Path(depth_file), genome_groups, roles))
        out.update(parse_flagstat(Path(flagstat_file)))
        dicts.append(out)

    df = pd.DataFrame(dicts)
    df.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()
