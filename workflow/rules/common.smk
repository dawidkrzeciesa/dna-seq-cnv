import glob
from os import path

import pandas as pd
import math

# read in sample sheet

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample_name": str, "target_bed": str},
        keep_default_na=False,  # keep empty entries empty strings, do not encode as NaN
    )
    .set_index(["sample_name"], drop=False)
    .sort_index()
)

### set global variables

TUMOR_SAMPLES = set(
    samples.loc[
        samples["alias"].str.startswith(config["alias_prefixes"]["tumor"]),
        "sample_name",
    ]
)

TUMOR_ALIASES = set(
    samples.loc[
        samples["alias"].str.startswith(config["alias_prefixes"]["tumor"]), "alias"
    ]
)

NORMAL_ALIASES = set(
    samples.loc[
        samples["alias"].str.startswith(config["alias_prefixes"]["normal"]), "alias"
    ]
)

### wildcard constraints


wildcard_constraints:
    sample="|".join(samples["sample_name"].drop_duplicates()),
    group="|".join(samples["group"].drop_duplicates()),
    tumor_alias="|".join(TUMOR_ALIASES),
    normal_alias="|".join(NORMAL_ALIASES),
    alias="|".join([*TUMOR_ALIASES, *NORMAL_ALIASES]),


### helper functions


def get_tumor_purity_setting(wildcards):
    if "tumor_purity" in samples.columns:
        purity = samples.loc[
            (samples["sample_name"] == wildcards.sample)
            & (samples["group"] == wildcards.group),
            "tumor_purity",
        ].squeeze()
        if not math.isnan(purity):
            return f"--purity {purity}"
    return ""


def get_chr_sex(wildcards):
    sex = get_sample_sex(wildcards)
    if sex == "male":
        return "-y"
    else:
        return ""


def is_tumor_sample_in_group(tumor_sample, group):
    tumor_group = samples.loc[
        (samples["sample_name"] == tumor_sample)
        & (samples["group"] == group)
        & (samples["alias"].str.startswith(config["alias_prefixes"]["tumor"])),
        :,
    ]
    return len(tumor_group) == 1


def get_group_sample_type(group, sample_type):
    sample_type_prefix = config["alias_prefixes"][sample_type]
    return samples.loc[
        (samples["group"] == group)
        & (samples["alias"].str.startswith(sample_type_prefix)),
        "sample_name",
    ].squeeze()


def get_normal_alias_of_group(group):
    normal_alias = samples.loc[
        (samples["group"] == group)
        & (samples["alias"].str.startswith(config["alias_prefixes"]["normal"])),
        "alias",
    ].squeeze()
    if len(normal_alias) == 0:
        normal_alias = samples.loc[
            samples["alias"].str.startswith(
                config["alias_prefixes"]["panel_of_normals"]
            ),
            "alias",
        ].squeeze()
    if type(normal_alias) != str and len(normal_alias) > 1:
        raise ValueError(
            f"Ambiguous normal sample for group '{group}'. Found more than one normal alias:\n {normal_alias}"
        )
    if len(normal_alias) == 0:
        raise ValueError(
            f"No normal alias available for group '{group}', but this is needed for purity estimation\n"
            "with theta2, when you do not specify a tumor_purity in the samples.tsv file."
        )
    else:
        return normal_alias


def get_sample_sex(wildcards):
    sex = samples.loc[
        (samples["sample_name"] == wildcards.sample)
        & (samples["group"] == wildcards.group),
        "sex",
    ].squeeze()
    return sex


def get_tumor_sample_group_aliases_combinations():
    combinations = []
    for tumor_sample in TUMOR_SAMPLES:
        groups = samples.loc[samples["sample_name"] == tumor_sample, "group"]
        if len(groups) == 0:
            raise ValueError(
                f"Tumor sample '{tumor_sample}' has no group assigned in the sample sheet."
            )
        else:
            for g in groups:
                if is_tumor_sample_in_group(tumor_sample, g):
                    tumor_alias = samples.loc[
                        (samples["group"] == g)
                        & (samples["sample_name"] == tumor_sample)
                        & (
                            samples["alias"].str.startswith(
                                config["alias_prefixes"]["tumor"]
                            )
                        ),
                        "alias",
                    ].squeeze()
                    normal_alias = get_normal_alias_of_group(g)
                    combinations.append(
                        f"{tumor_sample}.{g}.{tumor_alias}.{normal_alias}"
                    )
    return combinations


### input functions


def get_cnvkit_batch_input(wildcards, sample_type="tumor", ext="bam"):
    sample_name = get_group_sample_type(wildcards.group, sample_type)
    if (len(sample_name) == 0) and (sample_type == "normal"):
        # fall back to a `panel_of_normals` alias, if no matched normal sample is present
        sample_name = samples.loc[
            samples["alias"].str.startswith(
                config["alias_prefixes"]["panel_of_normals"]
            ),
            "sample_name",
        ].squeeze()
        if len(sample_name) == 0:
            # if not even a `panel_of_normal` sample is available, stick the tumor sample
            # here to be able to automatically trigger the rule with an empty `--normal`
            # specification via the params.normal lambda function
            sample_name = get_group_sample_type(wildcards.group, "tumor")
    if len(sample_name) == 0:
        raise AssertionError(
            f"Group {wildcards.group} does not have a sample with whose alias has the tumor prefix '{config['alias_prefixes']['tumor']}'"
        )
    bam_prefix = config.get("bam_prefix", "")
    return f"{bam_prefix}results/recal/{sample_name}.{ext}"


def get_cnvkit_batch_targets(wildcards):
    file = samples.loc[
        (samples["sample_name"] == wildcards.sample)
        & (samples["group"] == wildcards.group),
        "target_bed",
    ].squeeze()
    return file if file else []


def get_cnvkit_call_input(wildcards):
    # no purity specified for this sample
    if len(get_tumor_purity_setting(wildcards)) == 0:
        return (
            f"results/cnvkit_batch/{wildcards.sample}.{wildcards.group}.{wildcards.tumor_alias}.{wildcards.normal_alias}.purity_adjusted.cns",
        )
    else:
        return (
            f"results/cnvkit_batch/{wildcards.sample}.{wildcards.group}.{wildcards.tumor_alias}.{wildcards.normal_alias}.cns",
        )


def get_varlociraptor_present_bcf(wildcards):
    event = config["cnvkit"]["joint_event"]
    return f"results/calls/{wildcards.group}.SNV.{event}.fdr-controlled.bcf"
