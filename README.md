# Title in progress

This script is for user to use:
* mzML converted parquet file
* Spectrum index idx.json file
* Fragpipe searched psm.tsv file
  
to generate a npz file that can store a scipy sparse matrix. 

## Functions
```
# idx.json and parquet file are usually generated in the same parquet directory
parquet_dir = 'your/file/directory/species/parquet'
psm_dir = "your/file/directory/species/Peptide_search"

directory_processing(parquet_dir, psm_dir)
```



## Command Line Usage

