import math

import pandas as pd
import Aim1
import PSSMDictionaries as pd
import pickle
import operator


st_pssm = Aim1.st_pssm
y_pssm = Aim1.y_pssm
st_pssm_dict = pd.st_pssm_dict
y_pssm_dict = pd.y_pssm_dict


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
# k = [uniprot id,sequence] and v = all kinases and scores for that sequence
def substrate_scores(dict):
    substrate_scores = {}
    for key in dict.keys():
        current_seq = dict.get(key)
        temp_dict = st_score_dict(current_seq)
        if temp_dict is not None:
            uniprot_id = key[:key.index(',')]
            k = str(uniprot_id) + ',' + str(current_seq)
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


# create a dictionary to rank all the scores
# argument should be a dictionary where k = (uniprot id,sequence) and v = dictionary of kinases and scores
# key = (uniprot id, sequence) and value is a list of kinases from highest to lowest score
def ranked_dict(d):
    ranked_dict = {}
    for key, value in d.items():
        value = sort_dict(value)
        kinases = list(value.keys())
        kinase_rank_dict = {}
        i = 1
        while i < 304:
            for kinase in kinases:
                kinase_rank_dict[i] = kinase
                i += 1
        ranked_dict[key] = kinase_rank_dict
    return ranked_dict


# transform ranked dictionary into a dataframe
def ranked_df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    df = df.reset_index()
    df[['ACC_ID', 'SEQUENCE']] = df['index'].str.split(',', expand=True)
    df = df.drop(['index'], axis=1)

    col_one = df.pop('ACC_ID')
    col_two = df.pop('SEQUENCE')
    df.insert(0, 'ACC_ID', col_one)
    df.insert(1, 'SEQUENCE', col_two)

    return df



#print(y_score_dict('EALKFYTDPS'))




