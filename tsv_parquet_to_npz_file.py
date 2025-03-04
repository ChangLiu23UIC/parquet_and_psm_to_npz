import time
import pandas as pd
import os
import numpy as np
import scipy.sparse as sp
from pathlib import Path
import glob
import multiprocessing as mp
from itertools import repeat


# Constants
LOG_BASE_M = np.log(1.00001)
LOG200 = np.log(200)
MAX_DIM = int((np.log(2000) - np.log(200)) / np.log(1.00001)) + 1


def tsv_with_parquet_assemble(psm_file:str, parquet_file:str, json_file:str, result_path:str):
    """
    Uses the psm.tsv, parquet and idx.json file to output a npz file
    :param psm_file:
    :param parquet_file:
    :param json_file:
    :return:
    """
    start_func = time.time()
    # File name
    base_name = "".join(Path(parquet_file).name.split(".")[:-1])

    # Read files
    df_tsv = pd.read_csv(psm_file, sep='\t', usecols=["Spectrum", "Peptide", "Retention", "Charge", "Observed M/Z"])
    df_parquet = pd.read_parquet(parquet_file, columns=['mz', 'int'])
    df_idx = pd.read_json(json_file, orient="index")

    # Get stars index from idx.json for indptr
    arr_i = df_idx.to_numpy()
    indptr = np.append(arr_i[..., 0], arr_i[-1, 1]).astype(int)

    # Need Spectrum for merging (index joining) and  the change of file_name
    df_tsv['File'] = df_tsv['Spectrum']
    df_tsv['Spectrum'] = df_tsv["Spectrum"].str.split(".").str[-2].astype(np.int32)

    # Parquet columns To numpy arrays
    spec_arr = df_parquet["mz"].to_numpy()
    inten_arr = df_parquet["int"].to_numpy()

    # Normalize the intensities
    inten_arr = np.log2(inten_arr)
    inten_arr = np.concatenate(
        [(1 + ((x := inten_arr[i[0]:i[1]]) - (xmin := x.min())) / (x.max() - xmin) * 254).astype(np.uint8)
         for i in arr_i])
    # Indexing m/z
    spec_arr = ((np.log(spec_arr) - LOG200) / LOG_BASE_M).astype(np.int32)

    # Get the CSR matrix with intensity as value, spec_index as col_index and index starting site as indptr
    sparse_matrix = sp.csr_matrix((inten_arr, spec_arr, indptr), shape=((len(indptr) - 1), MAX_DIM))
    df_tsv["sparse_matrix"] = [sparse_matrix.getrow(i - 1) for i in df_tsv["Spectrum"]]

    # Convert to NPZ
    np_data = {col: df_tsv[col].to_numpy() for col in df_tsv.columns}
    np.savez(os.path.join(result_path, f"{base_name}.npz"), **np_data)

    print(f"{base_name} npz file generated at {result_path} for {time.time() - start_func} seconds")


def directory_processing(parquet_dir, psm_dir):
    """
    Parallel Processing of the function on the files
    :param parquet_dir:
    :param psm_dir:
    :return:
    """
    # Input Directories
    parquet_files = sorted(glob.glob(fr"{parquet_dir}/*.parquet"))
    json_files = sorted(glob.glob(fr"{parquet_dir}/*.idx.json"))
    psm_files = sorted(glob.glob(fr"{psm_dir}/*.tsv"))

    # Get Output Directory
    result_path = os.path.join(Path(parquet_dir).parent, "Merged_result")

    [print(i) for i in list(zip(psm_files, parquet_files, json_files, repeat(result_path)))]

    # Parallel processing the jobs
    with mp.Pool() as pool:
        pool.starmap(tsv_with_parquet_assemble, zip(psm_files, parquet_files, json_files, repeat(result_path)) , chunksize=1)


if __name__ == '__main__':
    psm_dir = r"G:\Files_for_Yuchi\Zebrafish\Peptide_search"
    parquet_dir = r'G:\Files_for_Yuchi\Zebrafish\parquet'

    start = time.time()
    directory_processing(parquet_dir, psm_dir)
    print(time.time() - start)