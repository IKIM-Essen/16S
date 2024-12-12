import os
import pandas as pd
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import move, copy2
from datetime import date, datetime
import sys


sys.stderr = open(snakemake.log[0], "w")

# Creating a metadata sample-sheet, that holds additional information concerning the samples.
# Metadata information is needed to be able to create plots and for metadata-specific analysis.

config = snakemake.config

DATA_PATH = str(config["data"])
IN_PATH = str(config["input"])
incoming_files = []
# Running through all files in the IN_PATH and checking them
# Adding the files to a list, when they meet the requirements
for f in listdir(IN_PATH):
    if (
        path.isfile(path.join(IN_PATH, f))
        and "ndetermined" not in f
        and ".fastq.gz" in f
        and os.stat(IN_PATH + f).st_size > 100
    ):
        incoming_files.append(f)
    else:
        print(f, "not used")
metadata = pd.read_csv(str(snakemake.input), header=0, delimiter=",")
date = metadata["run_date"].iloc[1]
# print(date)
# Adding a directory for the specific date as subdirectory to the DATA_PATH
if not os.path.isdir(DATA_PATH):
    mkdir(DATA_PATH)
DATA_PATH += date
if not os.path.isdir(DATA_PATH):
    mkdir(DATA_PATH)
# Checking if some files are already present in the DATA_PATH, those are not copied
data_files = [f for f in listdir(DATA_PATH) if path.isfile(path.join(DATA_PATH, f))]
files_not_to_copy = [f for f in data_files if f in incoming_files]
files_to_copy = [f for f in incoming_files if f not in data_files]
for file in files_to_copy:
    if file.endswith(".fastq.gz") and not "ndetermined" in file:
        copy2(IN_PATH + file, DATA_PATH)
# print(files_to_copy)

# Reading the sample names and the metadata from the file-names and the metadata.csv file, that needs to be provided
files = os.listdir(DATA_PATH)
# print(files)
sample_list = []
for name in files:
    sample = name.split("_")[0]
    sample_list.append(sample)
sample_list = list(set(sample_list))
metadata = pd.read_csv(str(snakemake.input), header=0, delimiter=",")
metadata.columns = metadata.columns.str.lower()
# Test parameter names for critical characters
column_headers = metadata.columns.tolist()
columns_with_characters = [
    col for col in column_headers if any(char in col for char in ["/", ".", "-"])
]
if columns_with_characters:
    print(
        "There are unprocessable characters in the metadata headers. Please make sure to eliminate dot, slash or hyphen and try again. \n Please have a look at the following columns"
    )
    print(columns_with_characters)
    raise Exception(
        "There are unprocessable characters in the metadata headers. Please make sure to eliminate dot, slash or hyphen and try again."
    )
# Samples, that are not mentioned in the metadata.csv are excluded from the sample-metadata-sheet
for name in sample_list:
    if name not in metadata["sample_name"].tolist():
        print("No metadata was found for " + name)
        sample_list.remove(name)
        # Also remove the fastq file from folder?

# Error if names don't fit
metadata_name_list = metadata["sample_name"].tolist()
metadata_name_list.remove("#q2:types")
if len(sample_list) == len(metadata_name_list):
    if set(sample_list) != set(metadata_name_list):
        print(
            "There seems to be a discrepancy between fastq file and metadata. Please check the spelling in the metadata.txt."
        )
        raise Exception(
            "There seems to be a discrepancy between fastq file and metadata. Please check the spelling in the metadata.txt."
        )

# Replacing empty metadata-columns with 0 and after that removing them from the file
metadata.fillna(0, inplace=True)
metadata.rename(columns={"sample_name": "sample-ID"}, inplace=True)
metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)

# Creating a sample_info df, that holds paths to the fastq files and the date.
# Common functions can read information from sample_info.txt instead from the data directories
data = [metadata["sample-ID"], metadata["run_date"]]
header = ["sample-ID", "run_date"]
sample_info = pd.concat(data, axis=1, keys=header)
if "SampleData[PairedEndSequencesWithQuality]" in str(snakemake.params.datatype):
    sample_info["path1"] = ""
    sample_info["path2"] = ""
elif "SampleData[SequencesWithQuality]" in str(snakemake.params.datatype):
    sample_info["path1"] = ""
sample_info.set_index("sample-ID", inplace=True)
sample_info.drop(labels=["#q2:types"], axis=0, inplace=True)
i = 0
while i < len(sample_info.index):
    for file in files_to_copy:
        if sample_info.index[i] in file and len(sample_info.index[i]) == len(
            file.split("_")[0]
        ):
            if "SampleData[PairedEndSequencesWithQuality]" in str(
                snakemake.params.datatype
            ):
                if "R1" in file:
                    sample_info.loc[sample_info.index[i], "path1"] = (
                        DATA_PATH + "/" + file
                    )
                elif "R2" in file:
                    sample_info.loc[sample_info.index[i], "path2"] = (
                        DATA_PATH + "/" + file
                    )
            elif "SampleData[SequencesWithQuality]" in str(snakemake.params.datatype):
                if "R1" in file:
                    sample_info.loc[sample_info.index[i], "path1"] = (
                        DATA_PATH + "/" + file
                    )
        else:
            continue
    i = i + 1
if not os.path.exists(str(snakemake.output.sample_info)):
    sample_info.to_csv(snakemake.output.sample_info, sep=",", mode="w")
