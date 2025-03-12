from re import findall as refindall
from collections import defaultdict

# the data of each amino acid and its composition
aminoacid = {
    'I': 'C6H13NO2',
    'L': 'C6H13NO2',
    'K': 'C6H14N2O2',
    'M': 'C5H11NO2S',
    'F': 'C9H11NO2',
    'T': 'C4H9NO3',
    'W': 'C11H12N2O2',
    'V': 'C5H11NO2',
    'R': 'C6H14N4O2',
    'H': 'C6H9N3O2',
    'A': 'C3H7NO2',
    'N': 'C4H8N2O3',
    'D': 'C4H7NO4',
    'C': 'C3H7NO2S',
    'E': 'C5H9NO4',
    'Q': 'C5H10N2O3',
    'G': 'C2H5NO2',
    'P': 'C5H9NO2',
    'S': 'C3H7NO3',
    'Y': 'C9H11NO3'
}

monoisotopic = {
    'S': 31.972072,
    'C': 12.0000,
    'H': 1.007825,
    'O': 15.994915,
    'N': 14.003074
}


def molecular_weight(molecule: str) -> float:
    """
    the function to calculate the molecular weight of the individual amino acid

    :param molecule: string
    :return: mass of the molecule
    """

    return sum(
        monoisotopic[atom] * int(num or '1')
        for atom, num in refindall(r'([A-Z][a-z]*)(\d*)', molecule)  # Use the regular expression
    )


def molecular_composition(molecule: str) -> dict:
    """
    It gives the molecular count for each atom within an amino acid
    :param molecule: peptide chain string
    :return: a dictionary of peptide in the peptide chain { peptide: 1(present)}
    """
    elements_counts = refindall('([A-Z][a-z]*)(\d*)', molecule)
    return {element: int(count) if count else 1 for element, count in elements_counts}


def peptide_composition(peptide: str) -> str:
    """
    This will count the composition of a peptide with all the atoms present
    :param peptide: peptide chain string
    :return:dictionary of the molecule { atom: count }
    """
    peptide_composition = defaultdict(int)
    for aa in peptide:
        mf = molecular_composition(aminoacid[aa])
        for element, count in mf.items():
            peptide_composition[element] += count

    # Remove the waters when they do peptide bonds
    peptide_composition["H"] -= 2 * (len(peptide) - 1)
    peptide_composition["O"] -= len(peptide) - 1

    return peptide_composition


def protein_index(protein: str) -> dict:
    """
    this function is to turn the protein into the amino acid dictionary to know the count

    :param protein: peptide chain string
    :return: a dictionary of peptide in the peptide chain { peptide: count }
    """
    protein_list = []
    protein_list[:0] = protein
    protein_peptides = {}
    # add a new key if the amino acid haven't showed up before, else, add one to the previous count
    for i in protein_list:
        if i in protein_peptides:
            protein_peptides[i] += 1
        else:
            protein_peptides[i] = 1
    return protein_peptides


#
def protein_weight(protein_dict: str) -> float:
    """
    this function will calculate the protein weight based on the number of each amino acid and its mass
    :param protein_dict: protein dict of amino acids
    :return: total mass of the protein
    """
    protein_final = protein_index(protein_dict)
    mass = 0
    for i in protein_final:
        mass += molecular_weight(aminoacid[i]) * protein_final[i]

    mass -= (len(protein_dict) - 1) * 18.010565
    return mass


if __name__ == '__main__':
    print(protein_weight("PEPTIDE"))
