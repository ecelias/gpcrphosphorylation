import pandas as pd
import Aim1_G
import EvaluatePSSM_G as ep
import math
from collections import OrderedDict

st_kinase_list = ep.st_kinase_list
y_kinase_list = ep.y_kinase_list
all_kinases = st_kinase_list + y_kinase_list


# create a variant sequence given a sequence
def var_seq(sequence, variant_position, phosphosite_position, variant_residue):
    altered_pos = int(variant_position) - int(phosphosite_position)
    if altered_pos == 0:
        var_seq = sequence[:5] + variant_residue + sequence[6:]
        return var_seq
    elif altered_pos != 0:
        altered_index = 5 + altered_pos
        var_seq = sequence[:altered_index] + variant_residue + sequence[altered_index + 1:]
        return var_seq


# function to merge two dictionaries
def merge_dict(dict1, dict2):
    merged_dict = dict(dict1.items() | dict2.items())
    return merged_dict


# function that accepts a sequence and returns a dictionary with k = kinase and v = log2(product) of all PSSMs
# dictionary will be for all kinases in yaffe paper
def st_score_dict_mut(sequence):
    st_score_dict_mut = {}
    sequence = str(sequence)
    for kinase in all_kinases:
        st_current_site = -5
        st_seq_index = 0
        st_product = 1
        while st_current_site < 5:
            if st_seq_index > len(sequence) - 1:
                break
            if kinase in y_kinase_list:
                st_product = 1
            elif kinase in st_kinase_list and st_current_site != 0:
                st_current_residue = sequence[st_seq_index]
                if st_current_residue == "_":
                    st_product *= 1
                elif ep.get_pssm_st(kinase, st_current_site, st_current_residue) is None:
                    st_product *= 1
                    print(f"Warning: Residue {st_current_residue} at {st_current_site} does not have PSSM score")
                else:
                    st_product *= ep.get_pssm_st(kinase, st_current_site, st_current_residue)
            st_seq_index += 1
            st_current_site += 1
        st_score_dict_mut[kinase + '_m'] = math.log2(st_product)
    return st_score_dict_mut


def y_score_dict_mut(sequence):
    y_score_dict_mut = {}
    sequence = str(sequence)
    for kinase in all_kinases:
        y_current_site = -5
        y_seq_index = 0
        y_product = 1
        while y_current_site < 5:
            if y_seq_index > len(sequence) - 1:
                break
            if kinase in st_kinase_list:
                y_product = 1
            elif kinase in y_kinase_list and y_current_site != 0:
                y_current_residue = sequence[y_seq_index]
                if y_current_residue == "_":
                    y_product *= 1
                elif ep.get_pssm_y(kinase, y_current_site, y_current_residue) is None:
                    y_product *= 1
                    print(f"Warning: Residue {y_current_residue} at {y_current_site} does not have PSSM score")
                else:
                    y_product *= ep.get_pssm_y(kinase, y_current_site, y_current_residue)
            y_seq_index += 1
            y_current_site += 1
        y_score_dict_mut[kinase + '_m'] = math.log2(y_product)
    return y_score_dict_mut


def phosphomimetic_score_dict_st(sequence):
    phosphomimetic_score_dict_st = {}
    for kinase in all_kinases:
        current_site = -5
        seq_index = 0
        product = 1
        while current_site < 5:
            if seq_index > len(sequence) - 1:
                break
            if kinase in y_kinase_list:
                product = 0
            elif current_site != 0:
                product = 0
            seq_index += 1
            current_site += 1
        phosphomimetic_score_dict_st[kinase + '_m'] = product
    return phosphomimetic_score_dict_st


def phosphomimetic_score_dict_y(sequence):
    phosphomimetic_score_dict_y = {}
    for kinase in all_kinases:
        current_site = -5
        seq_index = 0
        product = 1
        while current_site < 5:
            if seq_index > len(sequence) - 1:
                break
            if kinase in st_kinase_list:
                product = 0
            elif current_site != 0:
                product = 0
            seq_index += 1
            current_site += 1
        phosphomimetic_score_dict_y[kinase + '_m'] = product
    return phosphomimetic_score_dict_y


def null_score_dict_st(sequence):
    null_score_dict_st = {}
    for kinase in all_kinases:
        current_site = -5
        seq_index = 0
        product = 1
        while current_site < 5:
            if seq_index > len(sequence) - 1:
                break
            if kinase in y_kinase_list:
                product = 0
            elif current_site != 0:
                product = 0
            seq_index += 1
            current_site += 1
        null_score_dict_st[kinase + '_m'] = product
    return null_score_dict_st


def null_score_dict_y(sequence):
    null_score_dict_y = {}
    for kinase in all_kinases:
        current_site = -5
        seq_index = 0
        product = 1
        while current_site < 5:
            if seq_index > len(sequence) - 1:
                break
            if kinase in st_kinase_list:
                product = 0
            elif current_site != 0:
                product = 0
            seq_index += 1
            current_site += 1
        null_score_dict_y[kinase + '_m'] = product
    return null_score_dict_y


# for each mutated sequence, generate a product of PSSMs
# compare the total PSSM of WT and MT
# create a dictionary where k = [uniprot_id, mutated site] and v = score
def score_vars(df):
    var_substrate_scores = {}

    for index, row in df.iterrows():
        current_id = row['UNIPROT_ID']
        wt_seq = row['SEQUENCE']
        current_psite = row['PHOSPHOSITE']
        psite_residue = row['PS_RESIDUE']
        var_pos = row['VAR_POSITION']
        wt_res = row['WT_RESIDUE']
        mut_res = row['MUT_RESIDUE']
        mut_sequence = var_seq(wt_seq, var_pos, current_psite, mut_res)
        allele_count = row['ALLELE_COUNT']
        allele_num = row['ALLELE_NUM']
        allele_freq = row['ALLELE_FREQ']
        num_homo = row['NUM_HOMOZYG']

        wt_scores = {}
        mut_scores = {}
        if psite_residue.upper() in {'S', 'T'}:
            wt_scores = ep.st_score_dict(wt_seq)
        elif psite_residue.upper() == 'Y':
            wt_scores = ep.y_score_dict(wt_seq)

        if var_pos == current_psite and mut_res == 'Y':
            mut_scores = y_score_dict_mut(mut_sequence)
        elif var_pos == current_psite and mut_res in {'S', 'T'}:
            mut_scores = st_score_dict_mut(mut_sequence)

        elif psite_residue.upper() in {'S', 'T'} and var_pos != current_psite:
            mut_scores = st_score_dict_mut(mut_sequence)
        elif psite_residue.upper() == 'Y' and var_pos != current_psite:
            mut_scores = st_score_dict_mut(mut_sequence)

        elif var_pos == current_psite and mut_res in {'D', 'E'} and psite_residue.upper() in {'S', 'T'}:
            mut_scores = phosphomimetic_score_dict_st(mut_sequence)
        elif var_pos == current_psite and mut_res in {'D', 'E'} and psite_residue.upper() == 'Y':
            mut_scores = phosphomimetic_score_dict_y(mut_sequence)

        elif var_pos == current_psite and mut_res not in {'Y', 'S', 'D', 'E', 'T'} and psite_residue.upper() in {'S',
                                                                                                                 'T'}:
            mut_scores = null_score_dict_st(mut_sequence)
        elif var_pos == current_psite and mut_res not in {'Y', 'S', 'D', 'E', 'T'} and psite_residue.upper() in 'Y':
            mut_scores = null_score_dict_y(mut_sequence)

        scores = merge_dict(wt_scores, mut_scores)
        ordered_scores = OrderedDict(sorted(scores.items()))
        scores_sorted = dict(ordered_scores)

        if wt_scores is not None:
            k = str(current_id) + ',' + wt_seq + ',' + str(var_pos) + ',' + wt_res + ',' + mut_res + ',' + str(current_psite) + ',' + psite_residue + ',' + str(allele_count) + ',' + str(allele_num) + ',' + str(allele_freq) + ',' + str(num_homo)
            var_substrate_scores[k] = scores_sorted

    return var_substrate_scores


# transform the substrate_scores dictionary to a dataframe
def score_df(dictionary):
    df = pd.DataFrame.from_dict(dictionary, orient='index')
    df = df.reset_index(drop=True)
    df[['ACC_ID', 'WT_SEQUENCE', 'VARIANT_POSITION', 'WT_RES', 'MUT_RES', 'PHOSPHOSITE', 'PHOSPHOSITE_RESIDUE',
        'ALLELE_COUNT', 'ALLELE_NUM', 'ALLELE_FREQ', 'NUM_HOMOZYG']] = df['index'].str.split(',', expand=True)
    df = df.drop(['index'], axis=1)

    col_one = df.pop('ACC_ID')
    col_two = df.pop('WT_SEQUENCE')
    col_three = df.pop('VARIANT_POSITION')
    col_four = df.pop('WT_RES')
    col_five = df.pop('MUT_RES')
    col_six = df.pop('PHOSPHOSITE')
    col_seven = df.pop('PHOSPHOSITE_RESIDUE')
    col_eight = df.pop('ALLELE_COUNT')
    col_nine = df.pop('ALLELE_NUM')
    col_ten = df.pop('ALLELE_FREQ')
    col_eleven = df.pop('NUM_HOMOZYG')
    df.insert(0, 'ACC_ID', col_one)
    df.insert(1, 'WT_SEQUENCE', col_two)
    df.insert(2, 'VARIANT_POSITION', col_three)
    df.insert(3, 'WT_RES', col_four)
    df.insert(4, 'MUT_RES', col_five)
    df.insert(5, 'PHOSPHOSITE', col_six)
    df.insert(6, 'PHOSPHOSITE_RESIDUE', col_seven)
    df.insert(7, 'ALLELE_COUNT', col_eight)
    df.insert(8, 'ALLELE_NUM', col_nine)
    df.insert(9, 'ALLELE_FREQ', col_ten)
    df.insert(10, 'NUM_HOMOZYG', col_eleven)

    pd.options.display.float_format = '{:.3f}'.format
    return df


#gnomad_scores_dict = score_vars(Aim1_G.gnomad_phos_muts)
#gnomad_scores = score_df(gnomad_scores_dict)
#Aim1_G.to_csv(gnomad_scores, 'gnomad_kinase_scores')
