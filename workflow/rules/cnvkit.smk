rule cnvkit_access:
    input:
        fasta=get_reference,
    output:
        "results/access-mappable.bed",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/acces/access.log",
    params:
        access_param=config["cnvkit"]["access_param"],
    shell:
        "(cnvkit.py access {input} {params.access_param} -o {output}) 2>{log}"


rule cnvkit_batch:
    input:
        tumor=get_cnvkit_batch_input,
        normal=lambda w: get_cnvkit_batch_input(w, sample_type="normal"),
        bai_T=lambda w: get_cnvkit_batch_input(w, ext="bai"),
        bai_N=lambda w: get_cnvkit_batch_input(w, sample_type="normal", ext="bai"),
        fasta=get_reference,
        targets=get_cnvkit_batch_targets,
        access=rules.cnvkit_access.output,
    output:
        cns="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.cns",
        cnr="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.cnr",
        cnn="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.cnn",
    params:
        folder=lambda wc, output: path.dirname(output.cns),
        basename=lambda wc, input: path.splitext(path.basename(input.tumor))[0],
        normal=lambda wc, input: input.normal if (input.normal != input.tumor) else " ",
        batch=config["cnvkit"]["batch"],
        chr_sex=get_chr_sex,
        targets=lambda wc, input: f"--targets {input.targets}" if input.targets else "",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.log",
    threads: 64
    shell:
        "(cnvkit.py batch {input.tumor} "
        "  --normal {params.normal} "
        "  {params.targets} "
        "  --fasta {input.fasta} "
        "  --output-reference {output.cnn} "
        "  --access {input.access} "
        "  --output-dir {params.folder} "
        "  --diagram "
        "  --scatter "
        "  -p {threads} "
        "  {params.chr_sex} "
        "  {params.batch}; "
        " mv {params.folder}/{params.basename}.cns {output.cns}; "
        " mv {params.folder}/{params.basename}.cnr {output.cnr}; "
        ") 2>{log}"


rule genotype_snvs_to_vcf:
    input:
        bcf=get_varlociraptor_present_bcf,
    output:
        vcf="results/theta2/{sample}.{group}.genotypes.vcf",
    log:
        "logs/theta2/{sample}.{group}.genotypes.vcf.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor genotype < {input.bcf} | bcftools view -O v -o {output.vcf}"


rule create_allele_depth_annotation_per_sample:
    input:
        vcf="results/theta2/{sample}.{group}.genotypes.vcf",
    output:
        tsv="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv",
    log:
        "logs/theta2/{sample}.{group}.{alias}.allele_depths.log",
    conda:
        "../envs/vembrane.yaml"
    shell:
        "vembrane table "
        "  --header none "
        "  'CHROM,POS,REF,ALT,"
        '",".join([ str(int(FORMAT["DP"]["{wildcards.alias}"] - round(FORMAT["DP"]["{wildcards.alias}"] * FORMAT["AF"]["{wildcards.alias}"]))),'
        'str(int(round(FORMAT["DP"]["{wildcards.alias}"] * FORMAT["AF"]["{wildcards.alias}"]))) ]).replace(" ","")'
        "' "
        "{input.vcf} "
        ">{output.tsv} "
        "2>{log} "


rule bgzip_and_tabix_allele_depths:
    input:
        tsv="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv",
    output:
        gz="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv.gz",
        tbi="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv.gz.tbi",
    log:
        "logs/theta2/{sample}.{group}.{alias}.bgzip_and_tabix.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "( bgzip {input.tsv}; "
        "  tabix -b 2 -e 2 {output.gz};"
        ") 2>{log}"


rule annotate_tumor_allele_depth:
    input:
        vcf="results/theta2/{sample}.{group}.genotypes.vcf",
        tsv="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv.gz",
        tbi="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv.gz.tbi",
    output:
        vcf="results/theta2/{sample}.{group}.genotypes.{alias}_ad.vcf",
    log:
        "logs/theta2/{sample}.{group}.genotypes.{alias}_ad.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "( bcftools annotate "
        "    --annotations {input.tsv} "
        "    --samples {wildcards.alias} "
        "    --columns CHROM,POS,REF,ALT,FORMAT/AD "
        "     --header-line '##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">' "
        "    {input.vcf} "
        "    >{output.vcf} "
        ") 2>{log}"


use rule annotate_tumor_allele_depth as annotate_normal_allele_depth with:
    input:
        vcf="results/theta2/{sample}.{group}.genotypes.tumor_ad.vcf",
        tsv="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv.gz",
        tbi="results/theta2/{sample}.{group}.{alias}.allele_depths.tsv.gz.tbi",
    output:
        vcf="results/theta2/{sample}.{group}.genotypes.{tumor_alias}_ad.{alias}_ad.vcf",
    log:
        "logs/theta2/{sample}.{group}.genotypes.{tumor_alias}_ad.{alias}_ad.log",


rule cnvkit_to_theta2:
    input:
        cns="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.cns",
        vcf="results/theta2/{sample}.{group}.genotypes.{tumor_alias}_ad.{normal_alias}_ad.vcf",
    output:
        interval_count="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.interval_count",
        tumor_snps="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.tumor.snp_formatted.txt",
        normal_snps="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.normal.snp_formatted.txt",
    log:
        "logs/theta2/{sample}.{group}.{tumor_alias}.{normal_alias}.input.log",
    conda:
        "../envs/cnvkit.yaml"
    params:
        tumor_basename=lambda wc: f"{wc.sample}.{wc.group}.tumor.snp_formatted.txt",
        normal_basename=lambda wc: f"{wc.sample}.{wc.group}.normal.snp_formatted.txt",
    shell:
        "(cnvkit.py export theta "
        "  --vcf {input.vcf} "
        "  --sample-id {wildcards.tumor_alias} "
        "  --normal-id {wildcards.normal_alias} "
        "  --output {output.interval_count} "
        "  {input.cns}; "
        "  mv {params.tumor_basename} {output.tumor_snps}; "
        "  mv {params.normal_basename} {output.normal_snps} "
        ") 2>{log}"


rule theta2_purity_estimation:
    input:
        interval_count="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.interval_count",
        tumor_snps="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.tumor.snp_formatted.txt",
        normal_snps="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.normal.snp_formatted.txt",
    output:
        res="results/theta2/{sample}.{group}.{tumor_alias}.{normal_alias}/{sample}.{group}.{tumor_alias}.{normal_alias}/.BEST.results",
    log:
        "logs/theta2/{sample}.{group}.{tumor_alias}.{normal_alias}.BEST.log",
    conda:
        "../envs/theta2.yaml"
    params:
        out_dir=lambda wc, output: path.dirname(output.res),
    threads: 8
    shell:
        "(RunTHetA.py {input.interval_count} "
        "  --TUMOR_FILE {input.tumor_snps} "
        "  --NORMAL_FILE {input.normal_snps} "
        "  --DIR {params.out_dir} "
        "  --BAF "
        "  --NUM_PROCESSES {threads} "
        "  --FORCE "
        ") &>{log}"


rule theta2_to_cnvkit:
    input:
        res="results/theta2/{sample}.{group}.{tumor_alias}.{normal_alias}/{sample}.{group}.{tumor_alias}.{normal_alias}/.BEST.results",
    output:
        cns="results/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.purity_adjusted.cns",
    log:
        "logs/cnvkit_batch/{sample}.{group}.{tumor_alias}.{normal_alias}.purity_adjusted.log",
    conda:
        "../envs/cnvkit.yaml"
    params:
        out_dir=lambda wc, output: path.dirname(output.cns),
    shell:
        "(cnvkit.py import-theta "
        "  --output-dir {params.out_dir} "
        "  {output.cns} "
        "  {input.res} "
        ") 2>{log}"


rule cnvkit_call:
    input:
        cns=get_cnvkit_call_input,
    output:
        cns="results/cnvkit_call/{sample}.{group}.{tumor_alias}.{normal_alias}.cns",
    params:
        call_param=config["cnvkit"]["call_param"],
        tumor_purity=get_tumor_purity_setting,
        chr_sex=get_chr_sex,
        sample_sex=get_sample_sex,
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/call/{sample}.{group}.{tumor_alias}.{normal_alias}.cns.log",
    threads: 1
    shell:
        "(cnvkit.py call "
        "  {params.chr_sex} "
        "  -x {params.sample_sex} "
        "  --drop-low-coverage "
        "  -m clonal "
        "  {params.tumor_purity} "
        "  {input.cns} "
        "  -o {output.cns} "
        ") 2>{log}"


rule cnvkit_export_seg:
    input:
        cns=expand(
            "results/cnvkit_call/{sample_group_aliases}.cns",
            sample_group_aliases=get_tumor_sample_group_aliases_combinations(),
        ),
    output:
        "results/cnvkit/segments.seg",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/segments.log",
    threads: 1
    shell:
        "(cnvkit.py export seg {input.cns} -o {output}) 2>{log}"


rule gistic_get_input:
    input:
        "results/cnvkit/segments.seg",
    output:
        "results/cnvkit/gistic.input.seg",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/cnvkit/pandas.log",
    threads: 1
    script:
        "../scripts/gistic2_input.py"
