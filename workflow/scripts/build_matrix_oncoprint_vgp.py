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


# df = pd.concat([pd.read_csv(f,sep="\t").assign(SampleID=os.path.basename(f).split("_")[0]) for f in input_files])

df_amp = df
df_amp = df_amp[["gene", "cn", "SampleID"]]
df_amp["amp1"] = np.where(df_amp["cn"] == 3, "amp1", "")
df_amp["amp1"] = np.where(df_amp["cn"] == 4, "amp1", "")
df_amp["amp2"] = np.where(df_amp["cn"] >= 5, "amp2", "")
df_amp["amp3"] = np.where(df_amp["cn"] > 10, "amp3", "")
df_amp["del"] = np.where(df_amp["cn"] == 1, "LOSS", "")
df_amp["homdel"] = np.where(df_amp["cn"] < 1, "HOMDEL", "")
cols_amp = ["amp1", "amp2", "amp3", "del", "homdel"]
df_amp["amp"] = df_amp[cols_amp].apply(
    lambda row: ";".join(row.values.astype(str)), axis=1
)
df_amp = df_amp.drop(cols_amp, axis=1)

df_amp["temp"] = df_amp["SampleID"] + "_" + df_amp["gene"]
df_amp = df_amp.drop_duplicates(subset=["temp"], keep="last")

df_amp = df_amp.pivot(index="SampleID", columns="gene", values="amp")


with open(snakemake.output["matrix_vgp"], "w") as output_f:
    print(df_amp.to_csv(sep="\t", index=True), file=output_f)
