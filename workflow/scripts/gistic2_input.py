import pandas as pd

import sys

sys.stderr = open(snakemake.log[0], "w", buffering=1)

df = pd.read_csv(snakemake.input[0], sep="\t")

chromosome = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "X",
    "Y",
]

df = df[df["chrom"].isin(chromosome)]


df.to_csv(snakemake.output[0], sep="\t", index=False)
