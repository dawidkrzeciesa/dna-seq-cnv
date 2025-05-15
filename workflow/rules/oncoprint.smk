rule get_onco_kb:
    output:
        "resources/onco_kb/cancerGeneList.txt",
    log:
        "logs/onco_kb/get_cancerGeneList.log",
    params:
        onco_kb_download_link=config["onco_kb"],
    shell:
        "wget -O {output} {params.onco_kb_download_link} 2> {log}"


rule filter_oncogene:
    input:
        oncokb="resources/onco_kb/cancerGeneList.txt",
        cns="results/cnvkit_call/{sample}.{group}.{tumor_alias}.{normal_alias}.cns",
    output:
        cns_oncogene="results/cnvkit_call/filtered/oncogene/{sample}.{group}.{tumor_alias}.{normal_alias}_oncogene_only.cns",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_oncogene/{sample}.{group}.{tumor_alias}.{normal_alias}.log",
    threads: 1
    script:
        "../scripts/filter_oncogene.py"


rule filter_tumor_suppressor:
    input:
        oncokb="resources/onco_kb/cancerGeneList.txt",
        cns="results/cnvkit_call/{sample}.{group}.{tumor_alias}.{normal_alias}.cns",
    output:
        cns_tsg="results/cnvkit_call/filtered/tumor_supressor/{sample}.{group}.{tumor_alias}.{normal_alias}_tumor_supressor_only.cns",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_tumor_supressor/{sample}.{group}.{tumor_alias}.{normal_alias}.log",
    threads: 1
    script:
        "../scripts/filter_tumor_supressor.py"


rule filter_vgp:
    input:
        oncokb="resources/onco_kb/cancerGeneList.txt",
        cns="results/cnvkit_call/{sample}.{group}.{tumor_alias}.{normal_alias}.cns",
    output:
        cns_vgp="results/cnvkit_call/filtered/vgp/{sample}.{group}.{tumor_alias}.{normal_alias}_vgp.cns",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_vgp/{sample}.{group}.{tumor_alias}.{normal_alias}.log",
    threads: 1
    script:
        "../scripts/filter_vgp.py"


rule build_matrix_tumor_suppressor:
    input:
        expand(
            "results/cnvkit_call/filtered/tumor_supressor/{sample_group_aliases}_tumor_supressor_only.cns",
            sample_group_aliases=get_tumor_sample_group_aliases_combinations(),
        ),
    output:
        matrix_tsg="results/oncoprint/matrix/tumor_supressor_matrix.tsv",
    log:
        "logs/oncoprint/matrix/tumor_supressor_matrix.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_matrix_oncoprint_tumor_suppressor.py"


rule oncoprint_tumor_suppressors:
    input:
        matrix="results/oncoprint/matrix/tumor_supressor_matrix.tsv",
    output:
        report(
            "results/oncoprint/tumor_supressor_deletions_only.pdf",
            category="oncoprints",
            caption="../report/oncoprint_fusions.rst",
        ),
    log:
        "logs/oncoprint/tumor_supressor_deletions_only.log",
    conda:
        "../envs/oncoprint.yaml"
    params:
        title=config["oncoprint"]["title"],
    script:
        "../scripts/oncoprint_tumor_suppressor.R"


rule build_matrix_oncogene:
    input:
        expand(
            "results/cnvkit_call/filtered/oncogene/{sample_group_aliases}_oncogene_only.cns",
            sample_group_aliases=get_tumor_sample_group_aliases_combinations(),
        ),
    output:
        matrix_tsg="results/oncoprint/matrix/oncogene_matrix.tsv",
    log:
        "logs/oncoprint/matrix/oncogene_matrix.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_matrix_oncoprint_oncogene.py"


rule oncoprint_oncogene:
    input:
        matrix="results/oncoprint/matrix/oncogene_matrix.tsv",
    output:
        report(
            "results/oncoprint/oncogene_amplifications_only.pdf",
            category="oncoprints",
            caption="../report/oncoprint_fusions.rst",
        ),
    log:
        "logs/oncoprint/oncogene_amplifications_only.log",
    conda:
        "../envs/oncoprint.yaml"
    params:
        title=config["oncoprint"]["title"],
    script:
        "../scripts/oncoprint_oncogene.R"


rule build_matrix_vgp:
    input:
        expand(
            "results/cnvkit_call/filtered/vgp/{sample_group_aliases}_vgp.cns",
            sample_group_aliases=get_tumor_sample_group_aliases_combinations(),
        ),
    output:
        matrix_vgp="results/oncoprint/matrix/vgp_matrix.tsv",
    log:
        "logs/oncoprint/matrix/vgp_matrix.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_matrix_oncoprint_vgp.py"


rule oncoprint_vgp:
    input:
        matrix="results/oncoprint/matrix/vgp_matrix.tsv",
    output:
        report(
            "results/oncoprint/oncogene_vgp.pdf",
            category="oncoprints",
            caption="../report/oncoprint_fusions.rst",
        ),
    log:
        "logs/oncoprint/oncogene_vgp.log",
    conda:
        "../envs/oncoprint.yaml"
    params:
        title=config["oncoprint"]["title"],
    script:
        "../scripts/oncoprint_vgp.R"
