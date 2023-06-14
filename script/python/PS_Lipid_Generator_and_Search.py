import csv
import intervaltree
import mpmath
import numpy as np
import math
import json
from collections import defaultdict

def check_target_correctness(targets, ground_truth_examples):
    all_target_names = set()
    shorthand_to_formula = {}
    for ground_truth_example in ground_truth_examples:
        shorthand_to_formula[ground_truth_example['Species Shorthand']] = ground_truth_example['Formula'] 
    for target in targets:
        possible_names = [target["name"]] + target["isomers"]
        for name in possible_names:
            all_target_names.add(name)
            if name in shorthand_to_formula:
                if shorthand_to_formula[name] != target["formula"]:
                    raise Exception(target['name'] + " did not match ground truth")
    for ground_truth_example in ground_truth_examples:
        if ground_truth_example['Species Shorthand'] not in all_target_names:
            raise Exception(ground_truth_example['Species Shorthand'] + " not in generated lipids")
        

def process_feature_table(feature_table_csv_filepath, mass_tolerance_ppm=5):
    mz_interval_tree = intervaltree.IntervalTree()
    features = {}
    with open(feature_table_csv_filepath) as feature_table_fh:
        for feature in csv.DictReader(feature_table_fh, delimiter="\t"):
            mz = float(feature['mz'])
            delta_mz = mz / 1000000 * mass_tolerance_ppm
            mz_interval_tree.addi(mz - delta_mz, mz + delta_mz, feature)
            
            feature['rtime_lower'] = float(feature['rtime_left_base']) - 1 * (float(feature['rtime']) - float(feature['rtime_left_base']))
            feature['rtime_higher'] = 1 * (float(feature['rtime_right_base']) - float(feature['rtime'])) + float(feature['rtime_right_base'])

            features[feature['id_number']] = feature
    return features, mz_interval_tree

def read_target_list(target_list_csv):
    targets = []
    with open(target_list_csv) as target_list_fh:
        for target in csv.DictReader(target_list_fh):
            targets.append(target)
    return targets

def calculate_mass(element_dict):
    monoisotope_mass_table = {
        "C": 12.0,
        "H": 1.00782503223,
        "O": 15.99491461957,
        "P": 30.97376199842, 
        "N": 14.00307400443,
        "e": 0.000548579909065
    }

    element_masses = []
    for element, count in element_dict.items():
        if count > 0:
            element_masses.append(mpmath.fprod([count, monoisotope_mass_table[element]]))
        elif count == 0:
            pass
        else: 
            return None
    return mpmath.fsum(element_masses)

def add_formula_dicts(d1, d2):
    new_dict = defaultdict(int)
    for key, value in d1.items():
        new_dict[key] += value
    for key, value in d2.items():
        new_dict[key] += value
    return new_dict

def generate_formula(element_dict):
    formula = ''
    for element in ["C", "H", "N", "O", "P"]:
        count = element_dict[element]
        if count > 1:
            formula += element + str(count)
        elif count == 1:
            formula += element
        elif count == 0:
            pass
        else:
            return None
    return formula

# PS O-28:0	665.463172	1.007825032	664.455347	C34H68NO9P
# LPS O-28:1;O	665.463172	1.007825032	664.455347	C34H68NO9P
# LPS 28:0 665.46316977062 C34H68NO9P

def generate_LPSO(headgroup, max_mz=1000, min_mz=0):
    LPSO_modification = {"O": -2, "H": 4}
    LPSO_headgroup = add_formula_dicts(PS_headgroup, LPSO_modification)
    return generate_PS(LPSO_headgroup, name="LPS O-")

def generate_LPSO_O(headgroup, max_mz=1000, min_mz=0):
    LPSO_O_modification = {"O": -1, "H": 4}
    LPSO_O_headgroup = add_formula_dicts(PS_headgroup, LPSO_O_modification)
    new_lipids = []
    for lipid in generate_PS(LPSO_O_headgroup, name="LPS O-"):
        lipid["name"] += ";O"
        new_lipids.append(lipid)
    return new_lipids

def generate_PSO(headgroup, max_mz=1000, min_mz=0):
    PSO_modification = {"O": -1, "H": 2}
    PSO_headgroup = add_formula_dicts(PS_headgroup, PSO_modification)
    return generate_PS(PSO_headgroup, name="PS O-")

def generate_LPS(headgroup, max_mz=1000, min_mz=0):
    LPS_modification = {"O": -1, "H": 2}
    LPS_headgroup = add_formula_dicts(PS_headgroup, LPS_modification)
    return generate_PS(LPS_headgroup, name = "LPS ") 

def generate_PS(headgroup, max_mz=1000, min_mz=0, name="PS "):
    all_PS =[]
    CH2_mass = float(calculate_mass({"C": 1, "H": 2}))
    for num_CH2_added in range(0, math.ceil(max_mz // CH2_mass)):
        for double_bonds in range(0, num_CH2_added):
            modification_element_dict = {
                "C": num_CH2_added,
                "H": 2 * (num_CH2_added - double_bonds)
            }
            PS_dict = add_formula_dicts(headgroup, modification_element_dict)
            PS_mass = calculate_mass(PS_dict)
            PS_formula = generate_formula(PS_dict)
            if PS_mass and PS_formula and min_mz < PS_mass < max_mz:
                PS_name = name + str(num_CH2_added) + ":" + str(double_bonds)
                all_PS.append({
                    "name": PS_name,
                    "neutral_mass": PS_mass,
                    "formula": PS_formula,
                    "formula_dict": PS_dict
                })

    return all_PS

def generate_adducts(lipids, adduct_name, addduct_formula_diff):
    for lipid in lipids:
        new_formula_dict = add_formula_dicts(lipid["formula_dict"], addduct_formula_diff)
        lipid[adduct_name] = calculate_mass(new_formula_dict)
    return lipids

def prep_for_serialization(target):
    target['formula_dict'] = dict(target['formula_dict'])
    target['neutral_mass'] = float(target['neutral_mass'])
    target['[M-H+e]'] = float(target['[M-H+e]'])
    return target

def combine_isomers(targets):
    isomers = {}
    for target in targets:
        if target['formula'] not in isomers:
            target['isomers'] = []
            isomers[target['formula']] = target
        else:
            isomers[target['formula']]['isomers'].append(target['name'])
    return isomers.values()

def search_features(mz_interval_tree, targets):
    C13_mass = 13.00335483507
    C12_mass = 12.0
    C13_C12_mass_diff = C13_mass - C12_mass
    samples_to_include = [
        # MG CHECK THIS!!!!!!!
        "MT_20230308_006",
        "MT_20230308_008",
        "MT_20230308_010",
        "MT_20230308_012",
        "MT_20230308_014",
        "MT_20230308_016",
        "MT_20230308_018",
        "MT_20230308_020",
        "MT_20230308_022",
        "MT_20230308_024",
        "MT_20230308_026",
        "MT_20230308_028",
        "MT_20230308_030"
    ]

    for target in targets:
        #target["hits_isotopologue_chain"] = {}
        target["mz_only_hits"] = []
        for mz_only_match in [x.data for x in mz_interval_tree.at(target["[M-H+e]"])]:
            #print(mz_only_match["id_number"])
            target["mz_only_hits"].append(mz_only_match["id_number"])
            at_least_one_isopologue = True
            C13_count = 0
            isotopologues = [mz_only_match]
            #isotopologue_intensities.append([float(mz_only_match[sample_name]) for sample_name in samples_to_include])
            while at_least_one_isopologue:
                #print(C13_count)
                at_least_one_isopologue = False
                C13_count += 1
                isotopologue_mass_diff = mpmath.fprod([C13_C12_mass_diff, C13_count])
                isotopologue_mass = mpmath.fsum([target["[M-H+e]"], isotopologue_mass_diff])
                posssible_iso_matches = [x.data for x in mz_interval_tree.at(isotopologue_mass)]
                rt_filtered_possible_iso_matches = [iso for iso in posssible_iso_matches if float(mz_only_match['rtime_lower']) < float(iso['rtime']) < float(mz_only_match['rtime_higher'])]
                intensity_rt_filtered_possible_iso_matches = []
                for iso in rt_filtered_possible_iso_matches:
                    matching_intensity_counts = 0
                    previous_isotopologue_counts = 0
                    breaks_intensity_counts = 0
                    for sample_name in samples_to_include:
                        previous_isotopologue_intensity = float(isotopologues[-1][sample_name])
                        this_isotopologue_intensity = float(iso[sample_name])
                        if previous_isotopologue_intensity > 0:
                            previous_isotopologue_counts += 1
                            if this_isotopologue_intensity < previous_isotopologue_intensity:
                                matching_intensity_counts += 1
                        else:
                            if this_isotopologue_intensity > 0:
                                breaks_intensity_counts += 1
                    if previous_isotopologue_counts > 0 and matching_intensity_counts / previous_isotopologue_counts > 0.9 and breaks_intensity_counts == 0:
                        intensity_rt_filtered_possible_iso_matches.append(iso)
                preferred_isotopologue = None
                lowest_mass_error = np.inf
                for plausable_iso_match in intensity_rt_filtered_possible_iso_matches:
                    mass_error = abs(isotopologue_mass - float(plausable_iso_match['mz']))
                    if mass_error < lowest_mass_error:
                        preferred_isotopologue = plausable_iso_match
                if preferred_isotopologue:
                    isotopologues.append(preferred_isotopologue)
                    at_least_one_isopologue = True
                    if "hits_isotopologue_chain" not in target:
                        target["hits_isotopologue_chain"] = {}
                    if mz_only_match["id_number"] not in target["hits_isotopologue_chain"]:
                        target["hits_isotopologue_chain"][mz_only_match["id_number"]] = {0: [mz_only_match["id_number"]]}
                    if C13_count not in target["hits_isotopologue_chain"][mz_only_match["id_number"]]:
                        target["hits_isotopologue_chain"][mz_only_match["id_number"]][C13_count] = []
                    target["hits_isotopologue_chain"][mz_only_match["id_number"]][C13_count].append(preferred_isotopologue['id_number'])

def create_annotated_feature_table(all_targets, features):
    features_to_targets = {}
    for feature in features:
        features_to_targets[feature] = []
    for target in all_targets:
        for feature in target["mz_only_hits"]:
            assign_name = ';'.join([x+"_m+13C0" for x in [target["name"]] + target["isomers"]])
            features_to_targets[feature].append(assign_name)
        if "hits_isotopologue_chain" in target:
            for _, chain in target["hits_isotopologue_chain"].items():
                for isotopologue, feature_list in chain.items():
                    if isotopologue > 0:
                        for feature in feature_list:
                            features_to_targets[feature].append(';'.join([x+"_m+13C" + str(isotopologue) for x in [target["name"]] + target["isomers"]]))
    return features_to_targets

if __name__ == '__main__':
    import sys

    PS_headgroup = {
        "C": 6,
        "H": 10,
        "O": 10,
        "N": 1,
        "P": 1
    }

    features, mz_interval_tree = process_feature_table(sys.argv[1])
    min_mz = min([float(f['mz']) for f in features.values()])
    max_mz = max([float(f['mz']) for f in features.values()])
    all_targets = []
    all_targets.extend(generate_PS(PS_headgroup, min_mz = min_mz, max_mz = max_mz))
    all_targets.extend(generate_LPS(PS_headgroup, min_mz = min_mz, max_mz = max_mz))
    all_targets.extend(generate_PSO(PS_headgroup, min_mz = min_mz, max_mz = max_mz))
    all_targets.extend(generate_LPSO(PS_headgroup, min_mz = min_mz, max_mz = max_mz))
    all_targets.extend(generate_LPSO_O(PS_headgroup, min_mz = min_mz, max_mz = max_mz))
    all_targets = combine_isomers(all_targets)
    all_targets = generate_adducts(all_targets, "[M-H+e]", {"H": -1, "e": 1})

    search_features(mz_interval_tree, all_targets)
    with open("./theoretical_PS_search_results_MG_6_14_2023_v1.json", 'w') as out_fh:
        json.dump([prep_for_serialization(x) for x in all_targets], out_fh, indent=4)
    ground_truth = read_target_list(sys.argv[2])
    check_target_correctness(all_targets, ground_truth)

    with open("./features_to_PS_annotations_MG_6_14_2023_v1.json", "w") as out_fh:
        json.dump(create_annotated_feature_table(all_targets, features), out_fh, indent=4)
