import numpy as np
import pandas as pd
from scipy.sparse import coo_array


npz = np.load(r'G:\Files_for_Yuchi\Zebrafish\Merged_result\NP11-TDr-BA8-O1005_7-2.npz', allow_pickle= True)
df = pd.DataFrame({item: npz[item] for item in npz.files})

# npz2 = np.load(r'G:\Files_for_Yuchi\Fruit_Fly\Merged_result\Sample02.npz', allow_pickle= True)
# df2 = pd.DataFrame({item: npz2[item] for item in npz2.files})