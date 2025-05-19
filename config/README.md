Please set the path to the directory containing the BAM files and the reference genome in the config.yaml file.

In the `samples.tsv` file, use the same sample names as in the dna-seq-varlociraptor pipeline, which was used for generating the bam files.
If you are using data from targeted sequencing (whole exome sequencing or some other hybrid capture panel), you need to provide a `target_bed` file with the panel targets in the `samples.tsv`.
If you leave this column empty, the workflow will default to the `cnvkit.py batch --mode wgs`, so it will assume whole genome sequencing data.
