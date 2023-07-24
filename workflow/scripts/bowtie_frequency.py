import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import os
from pylab import savefig
import sys


sys.stderr = open(snakemake.log[0], "w")

path = os.path.dirname(str(snakemake.input[0]))
print(path)

inputs = []
for file in os.listdir(path):
    if file.endswith(".txt"):
        inputs.append(os.path.join(path, file))
print(type(inputs))
print(len(inputs))
# concatenate all txt files in a file called merged_file.txt
with open(str(snakemake.output), 'w') as outfile:
    i = 0
    while i < len(inputs):
        print(inputs[i])
        print(i)
        with open(inputs[i], encoding="utf-8", errors='ignore') as infile:
            path_name = inputs[i].split("_")[0]
            name = path_name.split("/")[-1] + ": " + infile.read()
            print(name)
            outfile.write(name)
            i += 1