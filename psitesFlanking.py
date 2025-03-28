import itertools
import Aim1
import pickle

final_mapped_psites = Aim1.readCSV("mapped_psites_with_flanking_info.csv")


def generate_phosphosite_permutations(up_id, psite, sequence, psite_res, flanking_psites, num_flanking_psites):

    result_dict = {}

    flanking_residues = list(sequence)

    if num_flanking_psites==0:
        key = f"{up_id}_{psite}_{psite_res}_permutation0"
        result_dict[key] = sequence
        return result_dict
    
    flanking_psite_idxs = flanking_psites.split(";")
    flanking_psite_idxs = [int(i) for i in flanking_psite_idxs]

    # Collect positions of T, Y, S residues in the flanking region (excluding center phosphosite)
    variable_positions = [
        i for i, r in enumerate(flanking_residues)
        if r.upper() in {'T', 'Y', 'S'} and i != 5 and i in flanking_psite_idxs
    ]

    # Generate all combinations of lowercase (phosphosite) or uppercase (not) for those positions
    permutations = []
    for bools in itertools.product([True, False], repeat=len(variable_positions)):
        mutated = flanking_residues.copy()
        for i, is_phospho in zip(variable_positions, bools):
            mutated[i] = mutated[i].lower() if is_phospho else mutated[i].upper()
        # Set center phosphosite lowercase
        mutated[5] = psite_res
        permutations.append(''.join(mutated))

    perm_idx = 1
    for perm in permutations:
        key = f"{up_id}_{psite}_{psite_res}_permutation{perm_idx}"
        result_dict[key] = perm
        perm_idx += 1

    # Create dictionary with all permutations
    key = f"{up_id}_{psite}_{psite_res}"

    return result_dict

def all_flanking_permutations():
    all_flanking_perms = {}

    for index, row in final_mapped_psites.iterrows():
        current_upid = row['UniProtID']
        current_psite = row['Position']
        current_seq = row['Sequence']
        current_psite_res = row['Phosphosite']
        current_flank_psites = row['FlankPSiteIndices']
        current_n_flank_psites = row['FlankPSiteCount']

        current_permutations = generate_phosphosite_permutations(current_upid, current_psite, current_seq, 
                                                                 current_psite_res, current_flank_psites, current_n_flank_psites)

        all_flanking_perms = {**all_flanking_perms, **current_permutations}
    
    return all_flanking_perms

all_flanking_permutations()

flanking_perms = all_flanking_permutations()
pickle.dump(flanking_perms, open("pickle_all_psites_with_flank_perms.p", 'wb'))
proteome_psites = pickle.load(open("pickle_all_psites_with_flank_perms.p", "rb"))
print(proteome_psites)
