#!/usr/bin/env python
# coding: utf-8

# In[1]:

import argparse
import pandas as pd
import numpy as np
import os
import subprocess
from sklearn import preprocessing
import sys


# In[2]:


parser = argparse.ArgumentParser()

parser.add_argument('--meta', type=str, required=True)
parser.add_argument('--clinical', type=str, required=True)
parser.add_argument('--expr', type=str, required=True)
parser.add_argument('--mapping', type=str, required=True)
parser.add_argument('--out_dir', type=str, required=True)

args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")


output_dir = args.out_dir

meta_file = args.meta
clinical_file = args.clinical
expr_file = args.expr
mapping_file = args.mapping
# In[3]:


# minimum number of control samples required for each cancer type
min_control_samples = 30

# fraction of patients where gene expression needs to be non zero, else the gene will be removed
min_fraction_expressed_patients = 0.8

# because of the log2(x+0.001) transform ot the tpms, a value of -9.9658 corresponds to a tpm of 0
zero_value = -9.9658


# Read meta file

# In[4]:


meta = pd.read_csv(meta_file, sep = "\t", index_col = 0)
print(f"Number of patients in meta: {len(meta.index)}")


#  Read clinical file

# In[5]:


clinical = pd.read_csv(clinical_file, sep = "\t", index_col = 0)
print(f"Number of patients in clinical: {len(clinical.index)}")


# In[6]:


subtypes = list(set(clinical["cancer type abbreviation"]))
subtypes.sort()
clinical = clinical.join(meta)
print(f"Number of patients in merged clinical: {len(clinical.index)}")
clinical.head()


# Read expression data

# In[7]:


expr = pd.read_csv( expr_file, sep = "\t", index_col = 0)
print(f"Number of patients with expression data: {len(expr.columns)}")


# Check which patients are in clinical and have expression data -> filter expression data accordingly

# In[8]:


shared_patients = list(set(clinical.index) & set(expr.columns))
print(f"Number of patients with expression and clinical data: {len(shared_patients)}")
expr = expr[shared_patients]
shared_clinical = clinical.loc[shared_patients]


# Keep only cancer types with at least 30 control samples

# In[9]:


to_keep = []
print("\nAvailable samples for each cancer type")
for st in subtypes:
    dat = shared_clinical[shared_clinical["cancer type abbreviation"] == st]
    control = dat.sample_type == "Solid Tissue Normal"
    if control.sum()>=min_control_samples:
        to_keep.append(st)
    print(f'{st}: cases {(control == 0).sum()}, controls {control.sum()}')


# In[10]:


print(f"\nCancer types with at least {min_control_samples} control samples:")

case_count = 0
control_count = 0
for st in to_keep:
    dat = shared_clinical[shared_clinical["cancer type abbreviation"] == st]
    control = dat.sample_type == "Solid Tissue Normal"
    print(f'{st}: cases {(control == 0).sum()}, controls {control.sum()}')
    case_count += (control == 0).sum()
    control_count += (control == 1).sum()
print(f'Total: cases {case_count}, controls {control_count}')


# Filter clinical, expression and shared clinical for patients in cancer types with enough control samples

# In[11]:


clinical = clinical[clinical["cancer type abbreviation"].isin(to_keep)]
shared_clinical = shared_clinical[shared_clinical["cancer type abbreviation"].isin(to_keep)]
expr = expr[shared_clinical.index]
print(f"\nNumber of patients available for selected cancer types: {len(clinical)}")
print(f"Number of patients with expression data available for selected cancer types: {len(shared_clinical)}")


# Read gene id to symbol mapping file, filter and remove entries with non unique symbols

# In[12]:


mapping = pd.read_csv(mapping_file, sep = "\t", index_col = 0)
print(f"\nEntries in the mapping file: {len(mapping.index)}")

duplicate_symbols = set(mapping.loc[mapping.gene.duplicated()]["gene"])
print(f"Number of non unique symbols in the mapping file: {len(duplicate_symbols)}")

mapping = mapping[~mapping.gene.isin(duplicate_symbols)]
print(f"Remaining entries: {len(mapping.index)}")
mapping = mapping.gene.to_dict()


# Map expression data to gene symbols

# In[13]:


expr = expr.loc[mapping.keys()]
expr = expr.rename(mapping)


# In[14]:


print(f"Shape of mapped expression data: {expr.shape}")


# Label control samples with 0 and case samples with 1

# In[15]:


shared_clinical["condition"] = (shared_clinical.sample_type != "Solid Tissue Normal")*1


# Filter for expressed genes and save data

# In[ ]:


print(f"\nFiltering for genes expressed in at least {min_fraction_expressed_patients*100}% of patients")
for st in to_keep:
    
    # select expression data of a specific cancer type
    dat = shared_clinical[shared_clinical["cancer type abbreviation"] == st]
    st_expr = expr[dat.index]
    
    control = dat.sample_type == "Solid Tissue Normal"
    
    expressed = np.sum(st_expr.T == zero_value) <= (st_expr.shape[1]) * (1-min_fraction_expressed_patients)
    st_expr = st_expr.loc[expressed]
    
    # scale based on mean and sd of control samples
    scaler = preprocessing.StandardScaler().fit(st_expr.loc[:,control].T)
    st_expr_scaled = scaler.transform(st_expr.T).T
    st_expr_scaled = pd.DataFrame(st_expr_scaled, columns=st_expr.columns, index=st_expr.index)
    
    
    print(f'{st}: cases {(control == 0).sum()}, controls {control.sum()}, genes {st_expr.shape[0]}')

    
    out_dir_st = os.path.join(output_dir,st)
    subprocess.call(f'mkdir -p {out_dir_st}', shell=True)
    
    dat.to_csv(os.path.join(out_dir_st, "tpm_meta.csv"))
    st_expr_scaled.to_csv(os.path.join(out_dir_st,"tpm.csv"))


# In[ ]:




