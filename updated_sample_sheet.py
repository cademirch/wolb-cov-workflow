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
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
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
        
    # If infection is provided as colon values, split them
    if isinstance(infection_value, str) and ':' in infection_value:
        return [inf.strip().lower() for inf in infection_value.split(':')]
    
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
        Dictionary of read files found
    """
    all_files = []
    
    # Search all directories for matching files
    for directory in read_dirs:
        if not directory.exists():
            logger.warning(f"Directory not found: {directory}")
            continue
        
        # Use more precise matching to avoid partial matches like JJ_MW_1 matching JJ_MW_100
        # Look for files that match the exact sample ID followed by underscore or end of string
        matching_files = []
        for f in directory.glob(f"*{sample_id}*fastq*"):
            # Check if it's an exact match by looking for patterns like:
            # - exact match: "JJ_MW_1.fastq.gz" 
            # - match followed by underscore: "JJ_MW_1_L001_R1.fastq.gz"
            # - match followed by underscore and S+digits: "JJ_MW_1_S1_L001_R1.fastq.gz"
            file_name = f.name
            if (f"_{sample_id}_" in f"_{file_name}" or 
                f"_{sample_id}." in f"_{file_name}" or
                f"/{sample_id}_" in f"/{f}" or
                file_name == f"{sample_id}.fastq.gz" or
                file_name.startswith(f"{sample_id}_")):
                matching_files.append(f)
        
        all_files.extend(matching_files)
    
    # Group by R1/R2 and ensure all paths are fully resolved (no .. or relative components)
    r1_files = sorted([str(Path(f).resolve()) for f in all_files if "_R1_" in str(f)])
    r2_files = sorted([str(Path(f).resolve()) for f in all_files if "_R2_" in str(f)])
    
    # Check if we found matching pairs
    if len(r1_files) != len(r2_files):
        logger.warning(f"Uneven number of R1 ({len(r1_files)}) and R2 ({len(r2_files)}) files for {sample_id}")
    
    if not r1_files and not r2_files:
        logger.warning(f"No read files found for sample {sample_id}")
    elif len(r1_files) < 1 or len(r2_files) < 1:
        logger.warning(f"Missing R1 or R2 files for sample {sample_id}: {len(r1_files)} R1, {len(r2_files)} R2")
    
    return {"read1": r1_files, "read2": r2_files}


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
            logger.error(f"Required column '{col}' not found in CSV. Available columns: {', '.join(df.columns)}")
            sys.exit(1)
    
    for col in recommended_columns:
        if col not in df.columns:
            logger.warning(f"Recommended column '{col}' not found in CSV. Some information may be missing.")
    
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
        if "Date Collected" in df.columns and (pd.isna(row["Date Collected"]) or not row["Date Collected"]):
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
            "read1_files": read_files["read1"],
            "read2_files": read_files["read2"],
        }
        
        # Add other metadata columns from original CSV
        for col in df.columns:
            if col != sample_id_col and col not in sample_entry:
                sample_entry[col] = row[col]
                
        sample_data.append(sample_entry)
    
    return pd.DataFrame(sample_data)


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Generate sample sheet for snakemake workflow")
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Path to the input CSV file with sample information"
    )
    parser.add_argument(
        "-d", "--directory",
        type=Path,
        action="append",
        required=True,
        dest="directories",
        help="Path to directory containing read files (can be specified multiple times)"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output JSON file path (defaults to stdout)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Process the directories to make them resolved absolute paths (removing any .. or relative components)
    resolved_dirs = [Path(directory).resolve() for directory in args.directories]
    
    # Process the data
    result_df = process_sample_data(args.input, resolved_dirs)
    logger.info(f"Processed {len(result_df)} samples")
    
    # Convert to JSON and output
    result_json = result_df.to_json(orient="records")
    
    if args.output:
        with open(args.output, "w") as f:
            f.write(result_json)
        logger.info(f"Output written to {args.output}")
    else:
        print(result_json)
        logger.info("Output written to stdout")


if __name__ == "__main__":
    main()