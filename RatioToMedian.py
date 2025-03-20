import PSSMDictionaries as pd
import EvaluatePSSM as ep
import Aim1
import statistics
import numpy as np

pssm = Aim1.pssm
pssm_dict = pd.pssm_dictionary
kinase_list = ep.kinase_list


# function to hold median values for each position
def flank_medians():
    median_dict = {}
    score_dict = {}
    temp = {}
    for index, row in pssm.iterrows():
        gene_name = row["GENE"]
        for column, item in pssm.items():
            if column != "GENE":
                residue_site = str(column)
                score = item[index]
                key = gene_name + ',' + residue_site
                score_dict[key] = score

    for key in score_dict.keys():
        gene = key[:key.index(',')]
        rs = key[key.index(',') + 1:]
        site = rs[:-1]
        key2 = gene + ',' + site
        score = score_dict.get(key)
        if key2 not in temp.keys():
            temp[key2] = [score]
        else:
            temp[key2].append(score)

    for key3 in temp.keys():
        current = temp.get(key3)
        median = statistics.median(current)
        median_dict[key3] = median
    return median_dict

# helper function to calculate fold change
def fold_change(score, median):
    l2fc = (score - median) / median
    return l2fc


# calculate fold change from median for flanking regions given a sequence
def flank_fc(seq):
    flank_score_dict = {}
    medians = flank_medians()
    seq = str(seq)
    for kinase in kinase_list:
        current_site = -5
        seq_index = 0
        fc = 0
        while current_site < 5:
            if seq_index > len(seq) - 1:
                break
            if current_site != 0:
                current_residue = seq[seq_index]
                current_pssm = ep.getPSSM(kinase, current_site, current_residue)
                key = str(kinase) + ',' + str(current_site)
                current_median = medians.get(key)
                if current_residue == "_":
                    fc += 0
                elif current_pssm is None:
                    fc += 0
                    print(f"Warning: Residue {current_residue} at {current_site} does not have PSSM score")
                else:
                    fc += fold_change(current_pssm, current_median)
            seq_index += 1
            current_site += 1
            flank_score_dict[kinase] = fc
    return flank_score_dict


test = flank_fc('QGGDFtRPNG')
print(test)

# calculate the S0/T0 values given a sequence
def favorability(sequence):
    S0_T0 = {}
    for kinase in kinase_list:
        current_site = -5
        seq_index = 0
        S = 0
        T = 0
        while current_site < 5:
            if current_site != 0:
                current_residue = sequence[seq_index]
                if current_residue == 'S':
                    S += ep.getPSSM(kinase, current_site, current_residue.lower())
                elif current_residue == 'T':
                    T += ep.getPSSM(kinase, current_site, current_residue.lower())
            seq_index += 1
            current_site += 1
        s_ctrl = 0.75 * S - 0.25 * T
        t_ctrl = 0.75 * T - 0.25 * S

        if not S == 0:
            S0 = s_ctrl / (max(s_ctrl, t_ctrl))
        else:
            S0 = 0
        if not T == 0:
            T0 = t_ctrl / (max(s_ctrl, t_ctrl))
        else:
            T0 = 0

        key_s = 'S0,' + str(kinase)
        key_t = 'T0,' + str(kinase)
        S0_T0[key_s] = S0
        S0_T0[key_t] = T0
    return S0_T0

# calculate median for S0/T0

