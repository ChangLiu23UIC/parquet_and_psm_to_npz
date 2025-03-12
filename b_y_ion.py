from protein_calculate import *
import pandas as pd


def cal_b_y_ion_mass(peptide):
    """
    Takes the peptide sequence and return a dataframe of b and y ion mass
    :param peptide: The peptide sequence
    :return: the dataframe of mass of b and y ions
    """
    sorted_list_of_mass = []
    protein_list = []
    protein_list[:0] = peptide.upper()
    b_and_y = []
    b_ions = []
    y_ions = []
    #split the peptides
    for i in range(1,len(protein_list)+1):
        # split the protein in different positions and group them into either b-ion or y-ion
        b_ion = protein_list[0:i]
        y_ion = protein_list[i - 1:]
        # the protein weight calculated previously contains extra H2O which needs to be removed while calculating the b or y ions
        # concatnate the list to make it suited for protein_weight function
        b_and_y.append(protein_weight(''.join(b_ion))+1 - 18.010565)
        b_and_y.append(protein_weight(''.join(y_ion))+19.010565 - 18.010565)
    # sort the list of mass
    sorted_list_of_mass = sorted(b_and_y)

    for i in range(0, len(sorted_list_of_mass)):
        if i % 2 == 0:
            b_ions.append(sorted_list_of_mass[i])
        else:
            y_ions.append(sorted_list_of_mass[i])

    indexing = list(range(1,len(peptide)+1))

    # Zip all the information into a table for pandas dataframe
    info_table =  zip(*[list(peptide.upper()), indexing, b_ions, y_ions[::-1], indexing[::-1]])

    df = pd.DataFrame(info_table, columns = ["Seq", "#", "B", "Y", "#(+1)"])

    b_fragments = [peptide[:i] for i in range(1, len(peptide))]
    y_fragments = [peptide[-i:] for i in range(1, len(peptide))]

    return df, b_fragments, y_fragments

if __name__ == '__main__':
    df, b_f, y_f = cal_b_y_ion_mass("PEPTIDE")
