1. The goal is to have a few mass spectra as input and predict the possible species, each given spectrum can refine the probability further, until we reach desired confidence (delta probability from first and second candidate species large enough), or when it reach stable state (probability won't change much). The data is already assembled into numpy data file (.npz) to facilitate the work and can be directly loaded.

   **Input**: few spectra

   ![image](https://www.broadinstitute.org/files/shared/proteomics/ms.jpg)

   **Output**: <img src="https://www.genesinspace.org/media/images/MOD_phylogenetic_tree.max-800x800.png" alt="image" style="zoom: 50%;" />

   

   

2. The npz file can be loaded as such:

   ```python
   import numpy as np
   import pandas as pd

   npz = np.load(r'G:\Files_for_Yuchi\Zebrafish\Merged_result\NP11-TDr-BA8-O1005_7-2.npz', allow_pickle= True)
   df = pd.DataFrame({item: npz[item] for item in npz.files})
   ```
Please notice that current npz files are separated by file and species

## Example 
The example dataframe structure read with pandas on the npz file. 
![image](https://github.com/user-attachments/assets/92e927cf-d65c-4987-822d-b39ea8a772a6)

   

3. Mass spectrum is an array of [(mass, intensity), (mass, intensity), ....], where mass monotonically increase, intensity, these are stored in the corresponding parquet file.

   **Example**: 

| MASS/CHARGE (M/Z) | Intensity |
| ----------------- | --------- |
| 103.1955          | 834.2     |
| 107.1497          | 958.5     |
| 149.4398          | 781.3     |
| 152.0573          | 2059.4    |
| 188.8318          | 934.9     |
| 195.4396          | 788.3     |
| 211.1402          | 1427.6    |
| 211.1496          | 4355.7    |
| 211.1622          | 1357.8    |
| 211.7111          | 1004.2    |
| 211.7215          | 8390.9    |
| 211.728           | 7958.2    |
| 211.7325          | 901.1     |
| 211.7377          | 2410.7    |
| 211.7425          | 1015.2    |
| 211.7469          | 2873.3    |
| 211.7513          | 2196.8    |
| 211.7556          | 2530.8    |
| 211.7599          | 1186      |
| 249.9105          | 951.2     |
| 302.2556          | 2459.3    |
| 355.5708          | 906.1     |
| 372.3964          | 942.6     |
| 374.3009          | 931.6     |
| 467.8401          | 1306.2    |
| 470.2251          | 960.9     |
| 495.835           | 10535.2   |
| 599.1929          | 955.1     |
| 660.5585          | 948       |
| 666.3334          | 1080.5    |
| 862.907           | 1144      |
| 923.229           | 1041.7    |
| 966.9485          | 983.7     |
| 984.098           | 1047.5    |
| 1182.922          | 1172.2    |
| 1256.552          | 1065.5    |

3. The mass/intensity array is stored as geometric progression sequence, with the understanding of:

   a. M/Z values are from 200 to 2000, anything above and below will be ignored

   b. M/Z values are very accurate, and contains most of the information, with an accuracy of +/- 5 parts per million (ppm), therefore, we constructed N = log(2000/200)/log(1+10ppm) = 230260 bins, and put the intensity to the corresponding bin, when two peaks are within 10ppm, we sum the intensity into the same bin

   c. The intensity is not as accurate and contains up to 15% uncertainty, therefore, for each spectrum, the intensity values are log transformed and normalized to int8, with value from 0 to 255

   d. After this transformation, we have obtained a sparse array of dimension (1, 230260), where the position index indicate the exact M/Z value, and the number inside is the log transformed and normalized intensity.

   e. This sparse matrix can be converted back to the mass/intensity array using the function below:

   ```python
   import numpy as np
   import scipy.sparse as sp

   LOG200 = np.log(200)
   LOG10PPM = np.log(1.00001)
   
   # mz_inten_csr is the sparse matrix <1 x 230260>
   mz_index = csr.indices
   inten = csr.values

   # convert the mz_index back to mz and combine with inten for an array of [(mass, intensity), (mass, intensity), ....]
   arr = np.stack((np.exp(mz_index * LOG10PPM + LOG200), inten), axis = -1)

   
   ```

   

4. raw parquet file can be read by the following function and the data structure is simple, the idx.json file contains the index of the start and end of each spectrum, for this project, only ms2 spectra are relevant:

   ```python
   import pandas as pd
   import b_y_ion
   

   df_parquet = pd.read_parquet(parquet_file)
   df_idx = pd.read_json(json_file, orient="index")

   
   # [add the ms2 display function here as well]
   ms2_df, b_frags, y_frags = cal_b_y_ion_mass("PEPTIDE")
   ```
   
   

5. In the data folder, parquet files are organized by different species, and each parquet file has a corresponding search results file (xxx_psm.tsv), which can be used as ground truth. This file is formatted as follow, we think only the charge, retention and precursor M/Z are useful for this project and these have been included in the npz file:

   | Spectrum                             | Spectrum File                                                | Peptide   | Modified Peptide            | Extended Peptide           | Prev AA | Next AA | Peptide Length | Charge |
   | ------------------------------------ | ------------------------------------------------------------ | --------- | --------------------------- | -------------------------- | ------- | ------- | -------------- | ------ |
   | NP11-TDr-BA8-O1005_7-2.00158.00158.2 | F:\Zebrafish_New\NP11_TDr_BA8_O1005_7_2\interact-NP11-TDr-BA8-O1005_7-2.pep.xml | SINSGGHK  |                             | ISAIVDGK.SINSGGHK.LGIGIEIE | K       | L       | 8              | 2      |
   | NP11-TDr-BA8-O1005_7-2.00202.00202.2 | F:\Zebrafish_New\NP11_TDr_BA8_O1005_7_2\interact-NP11-TDr-BA8-O1005_7-2.pep.xml | VGSAAQTR  |                             | NVGISVSR.VGSAAQTR.AMKQVAGT | R       | A       | 8              | 2      |
   | NP11-TDr-BA8-O1005_7-2.00218.00218.2 | F:\Zebrafish_New\NP11_TDr_BA8_O1005_7_2\interact-NP11-TDr-BA8-O1005_7-2.pep.xml | HPTDLDSSK | GYDPCNMK.HPTDLDSSK.IRGGMFDE | K                          | I       | 9       | 2              | 2      |
   | …                                    | …                                                            | …         | …                           | …                          | …       | …       | …              | …      |

   

6. We have already assembled the search results together with the parquet file into numpy data file (.npz), please see file xxx.py
