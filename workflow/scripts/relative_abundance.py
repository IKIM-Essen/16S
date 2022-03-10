import pandas as pd
import gzip
import shutil
import os
import zipfile
# Extracting the number of total features over all samples from the sample-table.
# Multiplying the number with the relative abundance filtering value to create a threshold for the 
# qiime filtering.

# Reading the sample-table, creating a zip file
file = str(snakemake.input)
name = os.path.splitext(file)[0]
shutil.copy(file, name + ".zip")
filename = name + ".zip"
# Extract zip files to folder
with zipfile.ZipFile(filename,"r") as zip_ref:
    name = filename.split("/")[-1]
    dir_name = os.path.dirname(str(snakemake.input))
    new_dir = dir_name + "/" + name
    zip_ref.extractall(os.path.splitext(new_dir)[0]+"/")
name = name.split(".")[0]
directory = os.path.dirname(str(snakemake.input)) + "/" + name
# Moving the folder inventory one folder up
b = 0
subdir = os.listdir(directory)
while b < len(subdir):
    orig_dir = directory + "/" + subdir[b]
    new_dir = directory
    for f in os.listdir(orig_dir):
        path = orig_dir + "/" + f
        shutil.move(path, new_dir)
    b = b + 1
# Read the specific csv holding the information, creating a dataframe, adding up all feature frequencies
datadir = directory +"/"+ "data/"
csv = datadir + "sample-frequency-detail.csv"
frequency = pd.read_csv(csv, header = 0, delimiter = ",")
number = frequency["0"].sum()
# Creating the abundance threshold and storing it in an output file
abundance = float(str(snakemake.params))
endnumber = number * abundance
endnumber = int(endnumber)
with open(str(snakemake.output), "w") as f:
    f.write(str(endnumber))


        