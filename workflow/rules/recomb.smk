rule find_chimeras:
    input:
        "results/bams/dedup/{sample}.bam",
    output:
        "results/recomb/chimeras/raw/{sample}.bam",
    script:
        "../scripts/find_recomb.py"


rule view_mappability:
    input:
        bam="results/recomb/chimeras/raw/{sample}.bam",
        mappability=get_mappability_bed
    output:
        filtered="results/recomb/chimeras/mappable/{e}/{sample}.bam",
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools index {input.bam}
        samtools view -h --region-file {input.mappability} {input.bam} | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -bh - |
        samtools sort -n - > {output.filtered}
        """

rule view_mappability2:
    input:
        bam="results/recomb/chimeras/raw/{sample}.bam",
        mappability=get_mappability_bed
    output:
        filtered="results/recomb/chimeras/mappable2/{e}/{sample}.bam",
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools index {input.bam}
        samtools view -h --region-file {input.mappability} {input.bam} | samtools view -bh - |
        samtools sort -n - > {output.filtered}
        """


rule get_rid_of_unpaired:
    input:
        bam="results/recomb/chimeras/mappable/{e}/{sample}.bam",
    output:
        bam="results/recomb/chimeras/final/{e}/{sample}.bam",
        bai="results/recomb/chimeras/final/{e}/{sample}.bam.bai",
        stats="results/recomb/chimeras/final/{e}/{sample}_algnstats.txt",
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools fixmate -u {input.bam} - | samtools view -b -e 'flag.paired' - |
        samtools sort -@ 5 -O BAM -o {output.bam} - ; samtools index -@ 5 {output.bam}
        samtools flagstat -O tsv {output.bam} > {output.stats}
        """


rule bamtobed:
    input:
        bam="results/recomb/chimeras/final/{e}/{sample}.bam",
        bai="results/recomb/chimeras/final/{e}/{sample}.bam.bai",
    output:
        bed="results/recomb/{e}/bambed/{sample}.bed",
    shell:
        "samtools sort -n {input.bam} | bedtools bamtobed -i stdin -bedpe > {output.bed}"