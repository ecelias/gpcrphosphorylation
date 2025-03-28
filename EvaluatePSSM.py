import math
import itertools
import pandas as pd
import Aim1
import pickle
import operator

from PhosphositeFeatures import all_mapped_psites
from PSSMDictionaries import st_pssm_dict
from PSSMDictionaries import y_pssm_dict


st_pssm = Aim1.st_pssm
y_pssm = Aim1.y_pssm


# create a list of all kinases in pssm
st_kinase_list = st_pssm['GENE'].tolist()
y_kinase_list = y_pssm['GENE'].tolist()
all_kinases = st_kinase_list + y_kinase_list


# helper function that accepts a kinase, site, and residue and returns the PSSM
def get_pssm_st(kinase, site, residue):
    key = str(kinase) + ',' + str(site) + ',' + str(residue)
    score = st_pssm_dict.get(key)
    if score is not None:
        score = float(score)
    return score


# pssm dict for tyrosine kinases
def get_pssm_y(kinase, site, residue):
    key = str(kinase) + ',' + str(site) + ',' + str(residue)
    score = y_pssm_dict.get(key)
    if score is not None:
        score = float(score)
    return score

# helper function to sort a dictionary of scores by score (value)
def sort_dict(d):
    sorted_dict = dict(sorted(d.items(), key=operator.itemgetter(1), reverse=True))
    return sorted_dict

# function that accepts a sequence (as a string) and returns a dictionary with k = kinase and v = log2(product) of all PSSMs
# dictionary will be for all kinases in yaffe paper
def st_score_dict(sequence):
    st_score_dict = {}
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
                elif get_pssm_st(kinase, st_current_site, st_current_residue) is None:
                    st_product *= 1
                    print(f"Warning: Residue {st_current_residue} at {st_current_site} does not have PSSM score")
                else:
                    st_product *= get_pssm_st(kinase, st_current_site, st_current_residue)
            st_seq_index += 1
            st_current_site += 1
        st_score_dict[kinase] = math.log2(st_product)
    return st_score_dict


def y_score_dict(sequence):
    y_score_dict = {}
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
                elif get_pssm_y(kinase, y_current_site, y_current_residue) is None:
                    y_product *= 1
                    print(f"Warning: Residue {y_current_residue} at {y_current_site} does not have PSSM score")
                else:
                    y_product *= get_pssm_y(kinase, y_current_site, y_current_residue)
            y_seq_index += 1
            y_current_site += 1
        y_score_dict[kinase] = math.log2(y_product)
    return y_score_dict


# create a dictionary where each substrate sequence has a score for each kinase
# k = uniprot id_sequence and v = all kinases and scores for that sequence
def substrate_scores():
    substrate_scores = {}
    for key in all_mapped_psites.keys():

        current_key = key.split(",")
        current_val = all_mapped_psites.get(key).split(",")

        current_seq = current_val[0]
        current_id = current_key[0]
        current_psite = current_key[1]
        current_psite_res = current_val[1]

        if current_psite_res in {"T", "S"}:
            temp_dict = st_score_dict(current_seq)
        elif current_psite_res == "Y":
            temp_dict = y_score_dict(current_seq)

        if temp_dict is not None:
            k = str(current_id) + '_' + str(current_seq) + "_" + current_psite_res
            substrate_scores[k] = temp_dict
    
    return substrate_scores


# transform the substrate_scores dictionary to a dataframe
def to_df(dictionary):
    df = pd.DataFrame.from_dict(dictionary, orient='index')
    df = df.reset_index()
    df[['ACC_ID', 'SEQUENCE']] = df['index'].str.split(',', expand=True)
    df = df.drop(['index'], axis=1)

    col_one = df.pop('ACC_ID')
    col_two = df.pop('SEQUENCE')
    df.insert(0, 'ACC_ID', col_one)
    df.insert(1, 'SEQUENCE', col_two)
    pd.options.display.float_format = '{:.3f}'.format
    return df




