import numpy as np
import pandas as pd
from protein_calculate import protein_weight


def cal_b_y_ion_mass(peptide):
    """
    Computes b and y ion masses for a given peptide sequence in a vectorized manner.
    :param peptide: The peptide sequence
    :return: A dataframe with b and y ion masses along with b and y ion fragments.
    """
    peptide = peptide.upper()
    n = len(peptide)

    # Generate b and y ion fragments
    b_fragments = [peptide[:i] for i in range(1, n + 1)]
    y_fragments = [peptide[-i:] for i in range(1, n + 1)]

    # Convert fragments into arrays
    b_ions = np.array([protein_weight(fragment) + 1.00784 - 18.010565 for fragment in b_fragments])
    y_ions = np.array([protein_weight(fragment) + 19.018405 - 18.010565 for fragment in y_fragments])

    # Construct dataframe
    df = pd.DataFrame({
        "Seq": list(peptide),
        "#": np.arange(1, n+1),
        "B": b_ions,
        "Y": y_ions[::-1],
        "#(+1)": np.arange(n , 0, -1)
    })

    return df, b_fragments[:-1], y_fragments[:-1]


if __name__ == '__main__':
    df, b_f, y_f = cal_b_y_ion_mass("PEPTIDE")
