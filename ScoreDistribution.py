import Aim1
import SubstratePhosphositeDictionary as spd
import EvaluatePSSM as ep
import EvaluateMut as em
import pickle
import pandas as pd

st_norm_kinase_list = ep.st_kinase_list
y_norm_kinase_list = ep.y_kinase_list

all_psites_dict = spd.all_mapped_psites
ht_psites_dict = spd.mapped_ht_psites


def dist_df(dictionary):
    df = pd.DataFrame.from_dict(dictionary, orient='index')
    df = df.reset_index()
    df[['ACC_ID', 'PHOSPHOSITE', 'SEQUENCE']] = df['index'].str.split(',', expand=True)
    df = df.drop(['index'], axis=1)

    col_one = df.pop('ACC_ID')
    col_two = df.pop('PHOSPHOSITE')
    col_three = df.pop('SEQUENCE')
    df.insert(0, 'ACC_ID', col_one)
    df.insert(1, 'SEQUENCE', col_three)
    df.insert(2, 'PHOSPHOSITE', col_two)
    pd.options.display.float_format = '{:.3f}'.format
    return df

# functions for ser/thr and tyr kinases
# for all mapped phosphosites get the current sequence and check what residue is at the phosphosite
# get the pssm score dictionary for each phosphosite sequence and duplicate for the mutant kinase
# merge wt and mut pssm dictionaries and return


def score_distribution_st(psites_dict):
    st_distribution_dict = {}
    for key in psites_dict.keys():
        print(key)
        current = psites_dict.get(key)
        current_seq = current[:current.index(',')]
        if current_seq[5] in {'S', 'T'}:
            all_pssm = ep.st_score_dict(current_seq)
            # pssm = ep.st_score_dict(current_seq)
            # pssm_mut = em.st_score_dict_mut(current_seq)
            # all_pssm = em.merge_dict(pssm, pssm_mut)
            k = key + ',' + current_seq
            st_distribution_dict[k] = all_pssm
    return st_distribution_dict


def score_distribution_y(psites_dict):
    y_distribution_dict = {}
    for key in psites_dict.keys():
        current = psites_dict.get(key)
        current_seq = current[:current.index(',')]
        if current_seq[5] in {'Y'}:
            print(current_seq)
            all_pssm = ep.y_score_dict(current_seq)
            #pssm = ep.y_score_dict(current_seq)
            #pssm_mut = em.y_score_dict_mut(current_seq)
            #all_pssm = em.merge_dict(pssm, pssm_mut)
            k = key + ',' + current_seq
            y_distribution_dict[k] = all_pssm
    return y_distribution_dict


def mean_dist_df(score_df, kinase_list):
    rows = ["mean", "std"]
    stats = pd.DataFrame(index=rows, columns=kinase_list)
    for kinase in kinase_list:
        mean = (score_df[kinase].astype(float)).mean()
        std = (score_df[kinase].astype(float)).std()
        stats.loc["mean", kinase] = mean
        stats.loc["std", kinase] = std
    return stats


# st_score_distribution = score_distribution_st(all_psites_dict)
# pickle.dump(st_score_distribution, open("pickle_st_score_distribution.p", 'wb'))
st_score_distribution_dict = pickle.load(open("pickle_st_score_distribution.p", "rb"))
st_score_dist = dist_df(st_score_distribution_dict)
#Aim1.to_csv(mean_dist_df(st_score_dist, st_norm_kinase_list), "st_stats_df")


# y_score_distribution = score_distribution_y(all_psites_dict)
# pickle.dump(y_score_distribution, open("pickle_y_score_distribution.p", 'wb'))
y_score_distribution_dict = pickle.load(open("pickle_y_score_distribution.p", "rb"))
y_score_dist = dist_df(y_score_distribution_dict)
#Aim1.to_csv(mean_dist_df(y_score_dist, y_norm_kinase_list), "y_stats_df")

all_psite_score_distribution = pickle.load(open("pickle_all_psite_score_distribution.p", "rb"))
all_psite_score_dist = dist_df(all_psite_score_distribution)

# ht_st_score_distribution = score_distribution_st(ht_psites_dict)
# pickle.dump(ht_st_score_distribution, open("pickle_ht_st_score_distribution.p", 'wb'))
ht_st_score_distribution_dict = pickle.load(open("pickle_ht_st_score_distribution.p", "rb"))
ht_st_score_dist = dist_df(ht_st_score_distribution_dict)
# Aim1.to_csv(ht_st_score_dist, "all_psp_ochoa_score_dist_st")
Aim1.to_csv(mean_dist_df(ht_st_score_dist, st_norm_kinase_list), "ht_st_stats_df")

# ht_y_score_distribution = score_distribution_y(ht_psites_dict)
# pickle.dump(ht_y_score_distribution, open("pickle_ht_y_score_distribution.p", 'wb'))
ht_y_score_distribution_dict = pickle.load(open("pickle_ht_y_score_distribution.p", "rb"))
ht_y_score_dist = dist_df(ht_y_score_distribution_dict)
# Aim1.to_csv(ht_y_score_dist, "all_psp_ochoa_score_dist_y")
Aim1.to_csv(mean_dist_df(ht_y_score_dist, y_norm_kinase_list), "ht_y_stats_df")

