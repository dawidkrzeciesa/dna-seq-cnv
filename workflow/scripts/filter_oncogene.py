import pandas as pd

import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

####### OncoKB #######

inport_columns = ["Hugo Symbol", "Is Oncogene", "Is Tumor Suppressor Gene"]

OncoKB = pd.read_csv(snakemake.input["oncokb"], sep="\t", usecols=inport_columns)

OncoKB = OncoKB.rename(
    columns={
        "Hugo Symbol": "gene",
        "Is Oncogene": "oncogene",
        "Is Tumor Suppressor Gene": "tsg",
    }
)

OncoKB = OncoKB[OncoKB.oncogene != "No"]

oncogene_lst = OncoKB["gene"].tolist()

####### filter cns #######
cns = pd.read_csv(snakemake.input["cns"], sep="\t")

cns["gene"] = cns.gene.str.split(",")
cns = cns.explode("gene").drop_duplicates()

cns_oncogene = cns[cns["gene"].isin(oncogene_lst)]

cns_oncogene.to_csv(snakemake.output["cns_oncogene"], sep="\t", index=False)
