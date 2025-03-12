import math
import itertools
import pandas as pd
from b_y_ion import *

def read_isotope_csv(filename: str) -> dict:
    """
    This will get a dictionary of abundances and mass for each of the isotope of the
    :param filename: str
    :return: dict
    """
    isotope_dict = {}
    isotopes = pd.read_csv("isotope.csv")
    for index, item in isotopes.iterrows():
        if item["Name"] in isotope_dict:
            isotope_dict[item["Name"]].append((item["Mass"], item["Abundances"]))
        else:
            isotope_dict[item["Name"]] = [(item["Mass"], item["Abundances"])]
    return isotope_dict


def get_combinations(isotope_num, atom_num):
    """
    Generate all combinations of isotopes that sum up to a certain atom number.

    Parameters:
    isotope_num (int): Number of isotopes/buckets to distribute atoms into.
    atom_num (int): Total number of atoms to distribute.

    Returns:
    list: List of all possible combinations of isotopes.
    """
    def generate_combinations(remaining_atoms, remaining_isotopes, current_combination):
        if remaining_isotopes == 0:
            if remaining_atoms == 0:
                return [current_combination]
            else:
                return []

        combinations = []
        for atoms in range(remaining_atoms + 1):
            combinations.extend(
                generate_combinations(
                    remaining_atoms - atoms,
                    remaining_isotopes - 1,
                    current_combination + [atoms]
                )
            )
        return combinations

    return generate_combinations(atom_num, isotope_num, [])


def multinomial_probability(n, outcomes):
    """
    Calculate the multinomial probability.

    :param n: Total number of trials (int)
    :param outcomes: A list of tuples, where each tuple contains the count and probability of each outcome.
    :return: The multinomial probability (float)
    """
    if sum([count for count, _ in outcomes]) != n:
        raise ValueError("The sum of the counts must equal n.")

    # Calculate the factorial of n
    factorial_n = math.factorial(n)

    # Calculate the product of the factorials of the counts and the power of probabilities
    denominator = 1
    probability_product = 1
    for count, probability in outcomes:
        if count < 0 or probability < 0 or probability > 1:
            raise ValueError("Counts and probabilities must be non-negative, and probabilities must not exceed 1.")
        denominator *= math.factorial(count)
        probability_product *= probability ** count

    # Calculate the multinomial probability
    multinomial_prob = factorial_n / denominator * probability_product
    return multinomial_prob


def sum_of_products(list1, list2):
    """
    Calculate the sum of products of corresponding elements from two lists.

    :param list1: First list of numbers
    :param list2: Second list of numbers
    :return: Sum of products of corresponding elements
    """
    if len(list1) != len(list2):
        raise ValueError("Both lists must have the same length.")

    # Calculate sum of products using a list comprehension and sum function
    result = sum(x * y for x, y in zip(list1, list2))
    return result


def prob_calc(atom, iso_num, total_num, iso_dict):
    """

    :param atom:
    :param iso_num:
    :param total_num:
    :param iso_dict:
    :return:
    """
    # Enumerate all the masses and possibliliteis for the isotope information
    prob_list = [mass[1] for mass in iso_dict[atom]]
    mass_list = [mass[0] for mass in iso_dict[atom]]
    diff_num = len(prob_list)-1
    prob_result_list = []
    mass_result_list = []
    for possible_isotopes in get_combinations(diff_num, iso_num):
        atom_distribution = [total_num-iso_num] + possible_isotopes
        if atom_distribution[0] <0:
            continue
        else:
            probability_result = multinomial_probability(total_num, list(zip(atom_distribution, prob_list)))
            prob_result_list.append(probability_result)

            mass = sum_of_products(atom_distribution, mass_list)
            mass_result_list.append(mass)

    return [prob_result_list, mass_result_list]


def prob_products(nested_list):
    """
    Compute all possible products from elements of each sublist in a nested list.

    :param nested_list: A list of lists, where each sublist contains numeric elements.
    :return: A list of products, one for each combination of elements from the sublists.
    """
    # Extract elements from sublists
    elements = [sublist for sublist in nested_list]

    # Compute the Cartesian product of these elements
    all_combinations = list(itertools.product(*elements))

    # Calculate the product for each combination
    products = [math.prod(combination) for combination in all_combinations]

    return products


def masses_sums(nested_list):
    """
    Compute all possible sums from elements of each sublist in a nested list.

    :param nested_list: A list of lists, where each sublist contains numeric elements.
    :return: A list of sums, one for each combination of elements from the sublists.
    """
    # Extract elements from sublists
    elements = [sublist for sublist in nested_list]

    # Compute the Cartesian product of these elements
    all_combinations = list(itertools.product(*elements))

    # Calculate the sum for each combination
    sums = [sum(combination) for combination in all_combinations]

    return sums

def isotope_calculator(peptide:str, iso_dict):
    """

    :param peptide:
    :param iso_dict:
    :return:
    """
    mass_key = []
    prob_val = []
    molecular_dict = peptide_composition(peptide)
    for i in range(5):
        comb_list = get_combinations(len(molecular_dict), i)
        atom_list = [atoms for atoms, count in molecular_dict.items()]
        for combs in comb_list:
            # will show how many of the atoms are there ex: [(4,C),(0,H)]
            features = list(zip(combs, atom_list))
            combs_dict = {}
            for iso_num, atom_name in features:
                probs_masses = prob_calc(atom_name, iso_num, molecular_dict[atom_name], iso_dict)
                combs_dict[atom_name] = probs_masses
            probs = 1
            mass = 0
            probs_prep = []
            masses_prep = []
            for atom_final, (probability, masses) in combs_dict.items():
                probs_prep += [probability]
                masses_prep += [masses]

            prob_res = prob_products(probs_prep)
            mass_res = masses_sums(masses_prep)


            for mass_name in mass_res:
                mass_key += [mass_name]
            for prob_name in prob_res:
                prob_val += [prob_name]

    results = dict(zip(mass_key, prob_val))
    filtered_dict = {key: value for key, value in results.items() if value > 0.01}
    sorted_dict = {key: filtered_dict[key] for key in sorted(filtered_dict)}
    final_result = average_and_sum_keys(sorted_dict, 0.01)

    return final_result


def average_and_sum_keys(data, threshold):
    sorted_keys = sorted(data)
    grouped_data = defaultdict(list)
    current_group = []

    for key in sorted_keys:
        if not current_group or key - current_group[-1] <= threshold:
            current_group.append(key)
        else:
            grouped_data[tuple(current_group)] = sum(data[k] for k in current_group)
            current_group = [key]

    if current_group:
        grouped_data[tuple(current_group)] = sum(data[k] for k in current_group)

    result = {sum(group) / len(group): value for group, value in grouped_data.items()}

    return result


if __name__ == '__main__':
    from b_y_ion import *

    isotope_dict = read_isotope_csv("isotope.csv")
    df, b_frag, y_frag = cal_b_y_ion_mass("SAMPLER")
