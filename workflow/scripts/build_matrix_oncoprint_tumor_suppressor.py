import pandas as pd
import numpy as np
import glob, os

input_files = str(snakemake.input)

input_files = input_files.split()

df = pd.concat(
    [
        pd.read_csv(f, sep="\t").assign(SampleID=os.path.basename(f).split("_")[0])
        for f in input_files
    ]
)


df_del = df[(df["cn"] < 2)]
df_del = df_del[["gene", "cn", "SampleID"]]
df_del["del"] = np.where(df_del["cn"] == 1, "LOSS", "HOMDEL")


df_del["temp"] = df_del["SampleID"] + "_" + df_del["gene"]

df_del = df_del.drop_duplicates(subset=["temp"], keep="last")

df_del = df_del.pivot(index="SampleID", columns="gene", values="del")


with open(snakemake.output["matrix_tsg"], "w") as output_f:
    print(df_del.to_csv(sep="\t", index=True), file=output_f)
