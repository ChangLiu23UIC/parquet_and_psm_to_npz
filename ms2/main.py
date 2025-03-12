from b_y_ion import *
import sys
import json
from flask import Flask

def main():
    peptide_sequence = sys.stdin.readline().strip()
    df = cal_b_y_ion_mass(peptide_sequence)

    print(df.to_json(orient="records"))


if __name__ == '__main__':
    main()
