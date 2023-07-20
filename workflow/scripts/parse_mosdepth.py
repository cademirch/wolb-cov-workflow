import pandas as pd
from pathlib import Path


def parse(file: Path, groups: dict[str : list[str]]) -> dict[str : str | str : float]:
    """
    Parses mosdepth summary file. Reports mean depth for groups described in groups arg.

    Args:
        file (Path): Path to mosdepth summary file
        groups (dict[str : list[str]]): Dict of coverage groups.
                            Key is group name, value is list of contigs

    Returns:
        dict[str : str | str : float]: {"sample_id":samp, "group1":0.0}
    """
    sample_id = file.name.replace(".mosdepth.summary.txt", "")
    inverse_groups = {}
    out = {}
    out["sample_id"] = sample_id
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
        out[group] = out[group] / len(contigs)
    return out


def create_output(
    df: pd.DataFrame, roles: dict[str : list[str]], outfile: Path
) -> None:
    """
    Writes csv of mean depths and titers.

    Args:
        df (pd.DataFrame): Dataframe with mean depth values
        roles (_type_): Dict describing host/symbiont
        outfile (Path): Path to output csv.
    """
    host = roles["host"]
    symbionts = roles["symbionts"]
    for symb in symbionts:
        df[f"{symb}_titer"] = df[symb] / df[host]
    df.to_csv(outfile, sep=",", index=False)


def main():
    groups = snakemake.params["groups"]  # noqa: F821
    roles = snakemake.params["roles"]  # noqa: F821
    infiles = snakemake.input  # noqa: F821
    outfile = Path(snakemake.output["csv"])  # noqa: F821
    depths = []

    for file in infiles:
        file = Path(file)
        depths.append(parse(file, groups))

    df = pd.DataFrame(depths)
    create_output(df, roles, outfile)


if __name__ == "__main__":
    main()
