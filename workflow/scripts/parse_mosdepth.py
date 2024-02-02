import pandas as pd
from pathlib import Path


def parse_flagstat(flagstat_file: Path, recomb: bool = False) -> dict:
    # Initialize variables to store mapped reads and duplicates
    mapped_reads = None
    duplicates = None
    total_reads = None

    try:
        # Open the flagstat file for reading
        with open(flagstat_file, "r") as file:
            # Read the lines from the file
            lines = file.readlines()

            # Loop through each line to find mapped reads and duplicates
            for line in lines:
                if "total" in line:
                    # Extract the number of total reads
                    total_reads = int(line.split()[0])
                elif "mapped" in line:
                    # Extract the number of mapped reads
                    mapped_reads = int(line.split()[0])
                elif "duplicates" in line:
                    # Extract the number of duplicates
                    duplicates = int(line.split()[0])

                # Break the loop if all values are found
                if all(
                    [x is not None for x in (mapped_reads, duplicates, total_reads)]
                ):
                    break

    except FileNotFoundError:
        print("File not found")

    # Return the extracted information as a dictionary
    if recomb:
        return {"recomb_reads": total_reads}
    return {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "duplicates": duplicates,
    }


def parse_depth(
    file: Path,
    groups: dict[str : list[str]],
    flagstat_files: list[Path],
    recomb_flagstat_files: list[Path],
) -> dict[str : str | str : float]:
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

    flagstat_file = Path([f for f in flagstat_files if sample_id in f][0])
    flagstats = parse_flagstat(flagstat_file)

    recomb_flagstat_file = Path([f for f in recomb_flagstat_files if sample_id in f][0])
    recomb_flagstats = parse_flagstat(recomb_flagstat_file, recomb=True)

    out.update(flagstats)
    out.update(recomb_flagstats)

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
    depths = snakemake.input.depths  # noqa: F821
    flagstat_files = snakemake.input.dedup_stats  # noqa: F821
    recomb_flagstat_files = snakemake.input.recomb_stats  # noqa: F821
    outfile = Path(snakemake.output["csv"])  # noqa: F821
    dicts = []
    for file in depths:
        file = Path(file)
        dicts.append(parse_depth(file, groups, flagstat_files, recomb_flagstat_files))

    df = pd.DataFrame(dicts)
    create_output(df, roles, outfile)


if __name__ == "__main__":
    main()
