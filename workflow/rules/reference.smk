rule combine_references:
    input:
        get_single_references
    output:
        "data/combined_genomes/{genome}.fa"
    shell:
        """
        cat {input} > {output}
        """


rule index_reference:
    input:
        ancient("data/combined_genomes/{genome}.fa")
    output:
        multiext("data/combined_genomes/{genome}.fa", ".fai", ".bwt", ".pac", ".ann", ".amb", ".sa"),
    log:
        "logs/index_reference/{genome}/log.txt",
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools faidx {input} &> {log}
        bwa index {input} &>> {log}
        """

rule bed_for_pileup:
    """
    Makes region bed of just wolbachia genomes for pileup
    """
    input:
        genome = ancient("data/combined_genomes/{genome}.fa.fai")
    output:
        bed = "data/bed/{genome}.bed"
    run:
        symbionts = config["roles"]["symbionts"]
        symbiont_chroms = []
        for k,v in config["coverage_groups"].items():
            if k in symbionts:
                symbiont_chroms.extend(v)
        print(symbiont_chroms)
        lines = []
        with open(input.genome, "r") as f:
            for line in f:
                line = line.strip().split()
                if line[0] in symbiont_chroms:
                    lines.append(line)
        print(lines)
        with open(output.bed, "w") as out:
            for l in lines:
                print(l[0], 0, l[1], sep="\t", file=out)




rule mappability_index:
    input:
        ancient("data/combined_genomes/{genome}.fa")
    output:
        index_dir=directory("data/mappability/{genome}/index"),
    log:
        "logs/mappability_index/{genome}/log.txt",
    conda:
        "../envs/env.yaml"
    shell:
        """
        genmap index -F {input} -I {output.index_dir} &> {log}
        """


rule mappability:
    input:
        index_dir=ancient("data/mappability/{genome}/index"),
    output:
        bedgraph=temp(
            "results/mappability/{genome}/{e}/{genome}.genmap.bedgraph"
        ),
        sorted_bedgraph="results/mappability/{genome}/{e}/sorted.bg",
        min_bed="results/mappability/{genome}/{e}/filtered.bed",
    log:
        "logs/mappability/{genome}/{e}/log.txt",
    params:
        k=config["mappability"]["kmer_size"],
        e="{e}",
        cutoff=config["mappability"]["score_cutoff"],
        out_dir=directory("results/mappability/{genome}/{e}"),
    conda:
        "../envs/env.yaml"
    shell:
        """
        genmap map -K {params.k} -E {params.e} -I {input.index_dir} -O {params.out_dir} \
            -bg -T {threads} -v &>> {log}
        sort -k1,1 -k2,2n {output.bedgraph} > {output.sorted_bedgraph} 2>> {log}
        awk 'BEGIN{{OFS="\\t";FS="\\t"}} {{ if($4>={params.cutoff}) print $1,$2,$3 }}' \
            {output.sorted_bedgraph} > {output.min_bed}
        """
