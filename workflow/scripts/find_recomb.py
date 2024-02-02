import pysam
from pathlib import Path

wmel_contig = "NC_002978.6"
wri_contig = "NC_012416.1"


def extract_chimeric_alignments(bam_file: pysam.AlignmentFile):
    out = []
    for read in bam_file.fetch(wmel_contig):
        if read.next_reference_name == wri_contig and not read.is_duplicate:
            out.append(read)
    for read in bam_file.fetch(wri_contig):
        if read.next_reference_name == wmel_contig and not read.is_duplicate:
            out.append(read)

    return out


def main():
    in_file = pysam.AlignmentFile(snakemake.input[0])  # noqa: F821
    algns = extract_chimeric_alignments(in_file)
    with pysam.AlignmentFile(
        filename=snakemake.output[0], mode="wb", template=in_file  # noqa: F821
    ) as out:
        for a in algns:
            out.write(a)


if __name__ == "__main__":
    main()
