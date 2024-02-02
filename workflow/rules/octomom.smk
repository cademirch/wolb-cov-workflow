rule mosdepth_octomom:
    input:
        bam="results/bams/filtered/{sample}.bam",
        bam_index="results/bams/filtered/{sample}.bam.bai",
    output:
        summary="results/octomom/{sample}.mosdepth.summary.txt",
        regions="results/octomom/{sample}.regions.bed.gz"
    log:
        "logs/mosdepth/octomom/{sample}.txt",
    params:
        prefix=lambda wc, output: output[0].replace(".mosdepth.summary.txt", ""),
    conda:
        "../envs/env.yaml"
    threads: resource_config["mosdepth"]["threads"]
    resources:
        mem_mb=resource_config["mosdepth"]["mem_mb"],
        runtime=resource_config["mosdepth"]["runtime"],
    shell:
        """
        mosdepth --by 1000 --chrom NC_002978.6 -F 1024 -n {params.prefix} {input.bam}
        """