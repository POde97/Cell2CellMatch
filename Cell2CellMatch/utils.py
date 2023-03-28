import os
import subprocess
import os
import torch
import random
import numpy as np 
from CellIDpy import ProteinCoding
import scanpy as sc

def runcmd(cmd, verbose = False, *args, **kwargs):

    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)
    pass
    
    
def DataDownloader():
  
  #os.makedirs("Cell2CellMatch")
  os.makedirs("data")
  os.makedirs("data/signature")
  os.makedirs("data/Glio")
  os.makedirs("data/Glio/Meta")
  os.makedirs("data/Glio/GBM")
  os.makedirs("data/Medullo")
  os.makedirs("data/Medullo/Meta")
  os.makedirs("data/Rabdo")
  os.makedirs("data/Rabdo/Meta")
  
  os.chdir("data/Medullo")
  #Medulloblastoma
  runcmd("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129730/suppl/GSE129730_RAW.tar")
  runcmd('tar -xf GSE129730_RAW.tar')
  runcmd('gunzip *gz') 
  os.chdir("..")
  #Glioblastoma
  os.chdir("Glio")
  runcmd("wget https://datahub-262-c54.p.genap.ca/GBM_paper_data/GBM_cellranger_matrix.tar.gz")
  runcmd('gunzip *gz')
  runcmd('tar -xf GBM_cellranger_matrix.tar -C GBM/')
  os.chdir("Meta")
  runcmd('wget --no-check-certificate https://docs.google.com/uc?export=download&id=1oi08Jjl89XoW2dejdh8wIZyM4lE1IN9e -O NeftelMeta.csv')
  runcmd('wget https://datahub-262-c54.p.genap.ca/GBM_paper_data/annotated_cancer_data.mat')
  os.chdir("..")
  runcmd('wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131928/suppl/GSE131928_RAW.tar')
  runcmd('tar -xf GSE131928_RAW.tar')
  runcmd('gunzip *gz') 
  runcmd('wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131928/suppl/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx')
  os.chdir("..")
  os.chdir("signature")
  runcmd('wget --no-check-certificate https://docs.google.com/uc?export=download&id=1GeJVEzBTxys9gTmV8caDaBffzNpkUlyd -O glioSign.csv')
  os.chdir("..")
  os.chdir("..")

  return "Data successfully downloaded"
  
  
def seed_everything(seed: int):
    r"""Sets the seed for generating random numbers in PyTorch, numpy and
    Python.

    Args:
        seed (int): The desired seed.
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def preproc(x):
    x = x[:,x.var_names.isin(ProteinCoding.HgProteinCoding())].copy()

    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
   

    return x


import pandas as pd
import csv
import os
import scipy.io
 

def load_ATAC(filename):
  mat = scipy.io.mmread(filename+ "matrix.mtx")
  peaks_path = filename+"peaks.bed"
  peaks=[]
  with open(peaks_path)as f:
    for line in f:
      L = line.strip().split()
      peaks.append(L[0].split(":")[0]+"_"+L[0].split(":")[1].split("-")[0]+"_"+L[0].split(":")[1].split("-")[1])
      #peaks.append(L[0])
  barcodes_path = filename+ "barcodes.tsv"
  barcodes=[]
  with open(barcodes_path)as f:
    for line in f:
      barcodes.append(line.strip().split()[0])
      
  adata = sc.AnnData(X = mat.T)
  adata.var_names = peaks
  adata.obs_names = barcodes
  adata.raw = adata.copy()


  return adata
  




