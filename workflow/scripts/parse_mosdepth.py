"""
Parses mosdepth summary output.
"""
from pathlib import Path
from collections import OrderedDict


def parse(file: Path, groups: OrderedDict[str : list[str]]) -> OrderedDict[str:float]:
    """
    Parses mosdepth summary file and returns dict
    where key is coverage group and value is avg coverage.
    """
    inverse_groups = OrderedDict()
    out = OrderedDict()
    for group, contigs in groups.items():
        out[group] = 0
        for contig in contigs:
            inverse_groups[contig] = group

    with open(file, "r") as f:
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
        out[group] = out[group] / len(contigs)
    return out


def main():
    groups = OrderedDict(snakemake.params["groups"])  # noqa: F821
    infiles = snakemake.input  # noqa: F821
    outfile = snakemake.output["csv"]  # noqa: F821

    with open(outfile, "w") as out:
        print("sample_id", groups.keys(), sep=",", file=out)
        for file in infiles:
            file = Path(file)
            sample_id = file.name.replace(".mosdepth.summary.txt", "")
            cov_dict = parse(file, groups)
            print(sample_id, cov_dict.values(), sep=",", file=out)


if __name__ == "__main__":
    main()
