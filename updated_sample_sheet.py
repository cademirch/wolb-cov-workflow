#!/usr/bin/env python3
"""
Makes sample sheet for snakemake workflow from CSV input.
Searches multiple read directories for matching files.
"""
import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Union

import pandas as pd

# Cell Line to host mapping
CELL_LINE_TO_HOST = {
    "s2": "dmel",
    "jw18": "dmel",
    "dsim6b": "dsim",
    "riv84": "dsim",
    "dsim": "dsim",
    "dmel": "dmel"
}

# Configure logging
def setup_logging(log_file=None, verbose=False):
    """Set up logging to both file and console"""
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Clear any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Create formatters
    verbose_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    simple_formatter = logging.Formatter("%(message)s")

    # Console handler (only warnings and errors to stderr)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(simple_formatter)
    root_logger.addHandler(console_handler)

    # File handler (all logs)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
        file_handler.setFormatter(verbose_formatter)
        root_logger.addHandler(file_handler)

    return root_logger


# Initialize logger (will be fully set up in main)
logger = logging.getLogger(__name__)


def process_infection(infection_value: str) -> Union[List[str], bool]:
    """
    Process infection value from CSV.
    
    Args:
        infection_value: Infection value from CSV
        
    Returns:
        List of infections or False if none
    """
    if pd.isna(infection_value) or not infection_value:
        return False

    # If infection is provided as comma-separated values, split them
    if isinstance(infection_value, str):
        if ":" in infection_value:
            return [inf.strip().lower() for inf in infection_value.split(":")]
        elif "+" in infection_value:
            return [inf.strip().lower() for inf in infection_value.split("+")]

    # Single infection
    return [infection_value.lower()]


def process_host(host_value: str) -> str:
    """
    Process host value from CSV using the global mapping.

    Args:
        host_value: Host value from CSV column "Cell Line"

    Returns:
        Mapped host value
    """
    if pd.isna(host_value) or not host_value:
        logger.warning("Missing host value in Cell Line column")
        return ""

    # Normalize the input by removing spaces and converting to lowercase
    normalized_host = host_value.lower().strip()

    # Check for direct match in the mapping
    if normalized_host in CELL_LINE_TO_HOST:
        return CELL_LINE_TO_HOST[normalized_host]

    # Check for partial matches
    for key in CELL_LINE_TO_HOST:
        if key in normalized_host:
            return CELL_LINE_TO_HOST[key]

    # If no match found, return the original value
    logger.warning(f"Couldn't map '{host_value}' to a known host. Using as is.")
    return host_value


def find_read_files(sample_id: str, read_dirs: List[Path]) -> Dict[str, List[str]]:
    """
    Find all read files for a sample across multiple directories.

    Args:
        sample_id: Sample identifier
        read_dirs: List of directories to search for reads

    Returns:
        Dictionary of read files found with keys 'read1_files' and 'read2_files'
    """
    all_files = []

    # Create alternate sample ID with colons converted to underscores
    # This is because sequencing facilities often convert special chars in filenames
    alt_sample_id = sample_id.replace(":", "_")

    logger.debug(
        f"Looking for read files for sample {sample_id} (and alternate {alt_sample_id})"
    )

    # Search all directories for matching files
    for directory in read_dirs:
        if not directory.exists():
            logger.warning(f"Directory not found: {directory}")
            continue

        # Use more precise matching to avoid partial matches like JJ_MW_1 matching JJ_MW_100
        # Look for files that match the exact sample ID followed by underscore or end of string
        matching_files = []

        # Try matching with the original sample ID
        for f in directory.glob(f"*{sample_id}*fastq*"):
            file_name = f.name
            if (
                f"_{sample_id}_" in f"_{file_name}"
                or f"_{sample_id}." in f"_{file_name}"
                or f"/{sample_id}_" in f"/{f}"
                or file_name == f"{sample_id}.fastq.gz"
                or file_name.startswith(f"{sample_id}_")
            ):
                matching_files.append(f)

        # Also try matching with the alternate sample ID (colons replaced with underscores)
        if sample_id != alt_sample_id:
            for f in directory.glob(f"*{alt_sample_id}*fastq*"):
                file_name = f.name
                if (
                    f"_{alt_sample_id}_" in f"_{file_name}"
                    or f"_{alt_sample_id}." in f"_{file_name}"
                    or f"/{alt_sample_id}_" in f"/{f}"
                    or file_name == f"{alt_sample_id}.fastq.gz"
                    or file_name.startswith(f"{alt_sample_id}_")
                ):
                    matching_files.append(f)

        all_files.extend(matching_files)

    # Group by R1/R2 and ensure all paths are fully resolved (no .. or relative components)
    r1_files = sorted([str(Path(f).resolve()) for f in all_files if "_R1_" in str(f)])
    r2_files = sorted([str(Path(f).resolve()) for f in all_files if "_R2_" in str(f)])

    # Check if we found matching pairs
    if len(r1_files) != len(r2_files):
        logger.warning(
            f"Uneven number of R1 ({len(r1_files)}) and R2 ({len(r2_files)}) files for {sample_id}"
        )

    if not r1_files and not r2_files:
        logger.warning(f"No read files found for sample {sample_id}")
    elif len(r1_files) < 1 or len(r2_files) < 1:
        logger.warning(
            f"Missing R1 or R2 files for sample {sample_id}: {len(r1_files)} R1, {len(r2_files)} R2"
        )

    return {"read1_files": r1_files, "read2_files": r2_files}


def process_sample_data(csv_file: Path, read_dirs: List[Path]) -> pd.DataFrame:
    """
    Process sample data from CSV and find associated read files.

    Args:
        csv_file: Path to input CSV file
        read_dirs: List of directories to search for reads

    Returns:
        DataFrame with processed sample data
    """
    # Read the CSV file
    try:
        df = pd.read_csv(csv_file)
        logger.info(f"Successfully read {csv_file} with {len(df)} entries")
    except Exception as e:
        logger.error(f"Failed to read CSV file {csv_file}: {e}")
        sys.exit(1)

    # Strip whitespace from string columns
    string_columns = df.select_dtypes(["object"])
    df[string_columns.columns] = string_columns.apply(lambda x: x.str.strip())

    # Ensure we have required columns
    required_columns = ["Sample ID"]
    recommended_columns = ["Cell Line", "Infection"]

    for col in required_columns:
        if col not in df.columns:
            logger.error(
                f"Required column '{col}' not found in CSV. Available columns: {', '.join(df.columns)}"
            )
            sys.exit(1)

    for col in recommended_columns:
        if col not in df.columns:
            logger.warning(
                f"Recommended column '{col}' not found in CSV. Some information may be missing."
            )

    sample_id_col = "Sample ID"

    # Process each sample
    sample_data = []

    for _, row in df.iterrows():
        sample_id = row[sample_id_col]

        # Skip empty sample IDs
        if pd.isna(sample_id) or not sample_id:
            logger.warning("Skipping row with empty Sample ID")
            continue

        # Skip rows with empty Date Collected
        if "Date Collected" in df.columns and (
            pd.isna(row["Date Collected"]) or not row["Date Collected"]
        ):
            logger.info(f"Skipping sample {sample_id} with empty Date Collected")
            continue

        # Find read files
        read_files = find_read_files(sample_id, read_dirs)

        # Extract host and infection from CSV columns
        host = process_host(row.get("Cell Line", ""))
        infection = process_infection(row.get("Infection", ""))

        # Create sample entry
        sample_entry = {
            "SampleID": sample_id,
            "host": host,
            "infection": infection,
            "read1_files": read_files["read1_files"],
            "read2_files": read_files["read2_files"],
        }

        # Add other metadata columns from original CSV
        for col in df.columns:
            if col != sample_id_col and col not in sample_entry:
                sample_entry[col] = row[col]

        sample_data.append(sample_entry)

    return pd.DataFrame(sample_data)


def main():
    """Main function"""
    # Track statistics for final summary
    stats = {
        "total_samples": 0,
        "samples_with_reads": 0,
        "samples_missing_reads": 0,
        "samples_missing_host": 0,
        "samples_missing_infection": 0,
        "samples_with_multiple_read_pairs": 0,
    }

    parser = argparse.ArgumentParser(
        description="Generate sample sheet for snakemake workflow"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        action="append",
        required=True,
        dest="input_files",
        help="Path to input CSV file with sample information (can be specified multiple times)",
    )
    parser.add_argument(
        "-d",
        "--directory",
        type=Path,
        action="append",
        required=True,
        dest="directories",
        help="Path to directory containing read files (can be specified multiple times)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="Output JSON file path (defaults to stdout)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "-l",
        "--log",
        type=Path,
        help="Log file path (defaults to samplesheet_prep_YYYYMMDD_HHMMSS.log)",
    )

    args = parser.parse_args()

    # Set up logging
    if not args.log:
        from datetime import datetime

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = Path(f"samplesheet_prep_{timestamp}.log")
    else:
        log_file = args.log

    setup_logging(log_file, args.verbose)
    logger.info(f"Starting sample sheet preparation, logging to {log_file}")

    # Process the directories to make them resolved absolute paths (removing any .. or relative components)
    resolved_dirs = [Path(directory).resolve() for directory in args.directories]
    logger.info(f"Searching for reads in {len(resolved_dirs)} directories")

    # Process each input CSV file and combine results
    all_sample_data = []
    samples_by_file = {}

    for input_file in args.input_files:
        logger.info(f"Processing input file: {input_file}")
        sample_df = process_sample_data(input_file, resolved_dirs)
        all_sample_data.append(sample_df)
        samples_by_file[str(input_file)] = len(sample_df)

    # Combine all dataframes
    if all_sample_data:
        result_df = pd.concat(all_sample_data, ignore_index=True)
        logger.info(
            f"Combined data from {len(args.input_files)} CSV files with total {len(result_df)} samples"
        )
        stats["total_samples"] = len(result_df)
    else:
        logger.error("No valid data found in any input files")
        sys.exit(1)

    # Calculate statistics
    for idx, row in result_df.iterrows():
        # Count samples with reads
        if row["read1_files"] and row["read2_files"]:
            stats["samples_with_reads"] += 1
            if len(row["read1_files"]) > 1 or len(row["read2_files"]) > 1:
                stats["samples_with_multiple_read_pairs"] += 1
                logger.info(
                    f"Sample {row['SampleID']} has multiple read pairs: {len(row['read1_files'])} R1, {len(row['read2_files'])} R2"
                )
        else:
            stats["samples_missing_reads"] += 1

        # Count samples missing host/infection
        if not row["host"]:
            stats["samples_missing_host"] += 1

        if not row["infection"] or row["infection"] is False:
            stats["samples_missing_infection"] += 1
    logger.info(f"Processed {len(result_df)} samples")

    # Check for duplicate sample IDs
    duplicate_samples = result_df["SampleID"].duplicated()
    if duplicate_samples.any():
        duplicates = result_df.loc[duplicate_samples, "SampleID"].tolist()
        logger.warning(
            f"Found duplicate sample IDs across input files: {', '.join(duplicates)}"
        )

    # Convert to JSON and output
    result_json = result_df.to_json(orient="records")
    
    if args.output:
        with open(args.output, "w") as f:
            f.write(result_json)
        logger.info(f"Output written to {args.output}")
    else:
        print(result_json)
        logger.info("Output written to stdout")

    # Print summary to stderr
    print(f"\nSample Sheet Summary:", file=sys.stderr)
    print(f"  Total samples processed: {stats['total_samples']}", file=sys.stderr)
    print(
        f"  Samples with reads: {stats['samples_with_reads']} ({stats['samples_with_multiple_read_pairs']} with multiple read pairs)",
        file=sys.stderr,
    )
    print(f"  Samples missing reads: {stats['samples_missing_reads']}", file=sys.stderr)
    print(f"  Samples missing host: {stats['samples_missing_host']}", file=sys.stderr)
    print(
        f"  Samples missing infection: {stats['samples_missing_infection']}",
        file=sys.stderr,
    )

    for input_file, count in samples_by_file.items():
        print(f"  {input_file}: {count} samples", file=sys.stderr)

    print(f"\nDetailed log available in: {log_file}", file=sys.stderr)


if __name__ == "__main__":
    main()