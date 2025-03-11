# Species recognition


This script is for user to use:
* mzML converted parquet file
* Spectrum index idx.json file
* Fragpipe searched psm.tsv file
  
to generate a npz file that can store a scipy sparse matrix. 

# Input
## Functions and inputs
### 1. Directory processing
```
# idx.json and parquet file are usually generated in the same parquet directory
parquet_dir = 'your/file/directory/species/parquet'
psm_dir = "your/file/directory/species/Peptide_search"

# Processes the entire species, generate npz files based on number of samples
directory_processing(parquet_dir, psm_dir)
```
### 2. Singe sample processing
```
parquet_file = 'your/file/directory/species/parquet/sample.parquet'
idx_file = 'your/file/directory/species/parquet/sample.idx.json'
psm_file = "your/file/directory/species/Peptide_search/sample.tsv"

# Process a single sample, generate one npz file
tsv_with_parquet_assemble(psm_file, parquet_file, idx_file) 
```

# Output
## Data Structure for NPZ file output 

**Spectrum** (From psm.tsv): The spectrum number

**Peptide** (From psm.tsv): The peptide sequence for the scan

**Charge** (From psm.tsv): The charge of the spectrum

**Retention** (From psm.tsv): The retention time of the spectrum 

**Observed M/Z** (From psm.tsv): The observed m/z value for the spectrum

**File** (From psm.tsv): The file name, the scan number and the charge of the spectrum

**Sparse_matrix_spec_index_intensity**
(Computed sparse matrix with psm.tsv and parquet):
The 1D sparse array with the indexed_m/z (m/z from 200 to 2000)  as the column and the normalized intensity(1 to 255) as the value.


## Read the saved NPZ file Output:
```
# example
import numpy as np
import pandas as pd

npz = np.load(r'G:\Files_for_Yuchi\Zebrafish\Merged_result\NP11-TDr-BA8-O1005_7-2.npz', allow_pickle= True)
df = pd.DataFrame({item: npz[item] for item in npz.files})
```


## Example 
The example dataframe structure read with pandas on the npz file. 
![image](https://github.com/user-attachments/assets/92e927cf-d65c-4987-822d-b39ea8a772a6)


