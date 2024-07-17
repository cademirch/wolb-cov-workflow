# Makes sample sheet for snakemake workflow from demux csv from duke
# Usage: python prepare_samplesheet.py <path to demux sheet> <path to reads dir> > <workflow_sheet>.json

from pathlib import Path
import argparse
import pandas as pd
import sys


def get_infection(samp: str) -> list[str] | bool:
    samp = samp.lower()
    choices = ["wmel", "wri", "wwil"]
    found = [inf for inf in choices if inf in samp]
    if found:
        return found
    print(f"infection not detected for {samp}, assuming none...", file=sys.stderr)
    return False


def get_host(samp: str) -> str:
    samp = samp.lower()
    dmel = ["s2", "jw18", "dmel"]
    dsim = ["dsim", "riv84", "dsin"]
    if any(d in samp for d in dmel):
        return "dmel"
    elif any(d in samp for d in dsim):
        return "dsim"

    print(f"couldnt figure out host for {samp}, assuming dmel", file=sys.stderr)
    return "dmel"


def get_read_names(
    sample_id: str, lane: int, sample_number: int, reads_dir: Path
) -> pd.Series:
    # {sample_id}_S{sample_number}_L{zfill(lane,3)}_R1_001.fastq.gz
    lane_num = str(lane).zfill(3)
    r1 = Path(
        reads_dir, f"{sample_id}_S{sample_number}_L{lane_num}_R1_001.fastq.gz"
    ).resolve()

    r2 = Path(
        reads_dir, f"{sample_id}_S{sample_number}_L{lane_num}_R2_001.fastq.gz"
    ).resolve()
    if not r1.exists():
        print(r1, "doesn't exist", file=sys.stderr)
        sys.exit(1)
    if not r2.exists():
        print(r2, "doesn't exist", file=sys.stderr)
        sys.exit(1)
    return pd.Series([str(r1), str(r2)], index=["read1", "read2"])


def get_read_names_no_sample_num(
    sample_id: str, lane: int, files: list[Path]
) -> pd.Series:

    found = [
        f.resolve() for f in files if sample_id in f.name and f"L00{lane}" in f.name
    ]
    if len(found) != 2:
        print(
            f"Expected 2 files but found {len(found)} for sample {sample_id} in lane {lane}: {found}",
            file=sys.stderr,
        )
        return None, None
    # Assuming the files are named such that read 1 and read 2 can be differentiated by R1 and R2
    read1 = [str(f) for f in found if "_R1_" in f.name][0]
    read2 = [str(f) for f in found if "_R2_" in f.name][0]
    return read1, read2


def main(file: Path, reads_dir: Path, ignore_sample_numbers: bool):
    # file must start with seq order id (4 digits)
    cols_to_keep = ["SampleID", "Lane", "# Reads", "Barcode"]
    df = pd.read_csv(file)
    valid_columns = [col for col in cols_to_keep if col in df.columns]
    df = df[valid_columns]
    string_columns = df.select_dtypes(["object"])
    df[string_columns.columns] = string_columns.apply(lambda x: x.str.strip())
    df = df.drop(df[df["SampleID"] == "Undetermined"].index)
    df["sequencing_order_id"] = file.stem[:4]

    if "# Reads" in df.columns:
        df = df[df["# Reads"] > 0]

    unique_sample_ids = df.drop_duplicates(subset="SampleID").reset_index(drop=True)

    unique_sample_ids["sample_number"] = range(1, len(unique_sample_ids) + 1)

    sample_number_mapping = unique_sample_ids.set_index("SampleID")[
        "sample_number"
    ].to_dict()

    df["sample_number"] = df["SampleID"].map(sample_number_mapping)
    if ignore_sample_numbers:
        files = list(reads_dir.glob("*.fastq.gz"))
        for index, row in df.iterrows():
            reads = get_read_names_no_sample_num(row["SampleID"], row["Lane"], files)
            df.at[index, "read1"] = reads[0]
            df.at[index, "read2"] = reads[1]
    else:
        df[["read1", "read2"]] = df.apply(
            lambda x: get_read_names(
                x["SampleID"], x["Lane"], x["sample_number"], reads_dir
            ),
            axis=1,
        )

    df["host"] = df["SampleID"].apply(get_host)
    df["infection"] = df["SampleID"].apply(get_infection)

    df.to_json(sys.stdout, orient="records")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "demux_csv_path",
        type=Path,
        help="Path to the demux CSV file",
    )
    parser.add_argument(
        "reads_dir_path",
        type=Path,
        help="Path to the directory where reads are stored",
    )
    parser.add_argument(
        "ignore_samp_num",
        type=bool,
        help="ignore sample nums for file matching",
        default=False,
    )

    args = parser.parse_args()

    main(args.demux_csv_path, args.reads_dir_path, args.ignore_samp_num)
