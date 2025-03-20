import pandas as pd
import EvaluateMut as em
from scipy.stats import norm
import Aim1
import pickle
import EvaluatePSSM as ep
import time

norm_kinase_list = em.all_kinases
st_norm_kinase_list = ep.st_kinase_list
y_norm_kinase_list = ep.y_kinase_list

sd_y = pd.read_csv('y_stats_df.csv', index_col=0)
sd_st = pd.read_csv('st_stats_df.csv', index_col=0)

def get_mut_kinases(kinases_list):
    mut_list = []
    for kinase in kinases_list:
        current_mut = kinase + '_m'
        mut_list.append(current_mut)
    mut_wt_kinase_list = kinases_list + mut_list
    mut_wt_kinase_list.sort()
    return mut_wt_kinase_list


kinase_list = get_mut_kinases(norm_kinase_list)
st_kinase_list = get_mut_kinases(st_norm_kinase_list)
y_kinase_list = get_mut_kinases(y_norm_kinase_list)


# function to sort a dictionary
def sorted_dict(dictionary):
    sorted_dict_by_values = sorted(dictionary.items(), key=lambda x: x[1], reverse=True)
    back_to_dict = dict(sorted_dict_by_values)
    return back_to_dict


# function to accept a sequence and kinase and return the percentile compared to
def percent_rank_dict(df):
    start_time = time.time()
    percent_dict = {}
    for index, row in df.iterrows():
        kinase_dict = {}
        print(index)
        for kinase in norm_kinase_list:
            mut_kinase = kinase + "_m"
            if kinase in st_norm_kinase_list:
                mean = sd_st.loc["mean", kinase]
                std = sd_st.loc["std", kinase]

                score = row[kinase]
                mut_score = row[mut_kinase]

                z_score = (score - mean) / std
                mut_z_score = (mut_score - mean) / std

                percent = norm.cdf(z_score) * 100
                mut_percent = norm.cdf(mut_z_score) * 100

                kinase_dict[kinase] = percent
                kinase_dict[mut_kinase] = mut_percent
            elif kinase in y_norm_kinase_list:
                mean = sd_y.loc["mean", kinase]
                std = sd_y.loc["std", kinase]

                score = row[kinase]
                mut_score = row[mut_kinase]

                z_score = (score - mean) / std
                mut_z_score = (mut_score - mean) / std

                percent = norm.cdf(z_score) * 100
                mut_percent = norm.cdf(mut_z_score) * 100

                kinase_dict[kinase] = percent
                kinase_dict[mut_kinase] = mut_percent

        key = str(row['ACC_ID']) + ',' + str(row['WT_SEQUENCE']) + ',' + str(row['VARIANT_POSITION']) + ',' + str(row['WT_RES']) + ',' + str(row['MUT_RES']) + ',' + str(row['PHOSPHOSITE']) + ',' + str(row['PHOSPHOSITE_RESIDUE'])
        percent_dict[key] = kinase_dict
        end = time.time()
        elapsed = end - start_time
        print(elapsed)
    return percent_dict


# transform percent rank dictionary into df
def percent_df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    df = df.reset_index()
    df[['ACC_ID', 'WT_SEQUENCE', 'VARIANT_POSITION', 'WT_RES', 'MUT_RES', 'PHOSPHOSITE', 'PHOSPHOSITE_RESIDUE']] = df['index'].str.split(',', expand=True)
    df = df.drop(['index'], axis=1)

    col_one = df.pop('ACC_ID')
    col_two = df.pop('WT_SEQUENCE')
    col_three = df.pop('VARIANT_POSITION')
    col_four = df.pop('WT_RES')
    col_five = df.pop('MUT_RES')
    col_six = df.pop('PHOSPHOSITE')
    col_seven = df.pop('PHOSPHOSITE_RESIDUE')
    df.insert(0, 'ACC_ID', col_one)
    df.insert(1, 'WT_SEQUENCE', col_two)
    df.insert(2, 'VARIANT_POSITION', col_three)
    df.insert(3, 'WT_RES', col_four)
    df.insert(4, 'MUT_RES', col_five)
    df.insert(5, 'PHOSPHOSITE', col_six)
    df.insert(6, 'PHOSPHOSITE_RESIDUE', col_seven)

    pd.options.display.float_format = '{:.3f}'.format

    return df


#yaffe = ep.yaffe_saved_scores
#yaffe_df = ep.to_df(yaffe)


# mut_scores = em.test2
# test3 = percent_rank_dict(mut_scores)
# test4 = percent_df(test3)
# Aim1.to_csv(test4, 'pr.test')
# print(test4)

#test_df = percent_df(test)
#pickle.dump(test_df, open("percentile_rank.p", "wb"))
#test11 = pickle.load(open("percentile_rank.p", "rb"))
#print(test11)

# scores = Aim1.cosmic_scores
# cosmic_rank_dict = percent_rank_dict(scores)
# cosmic_percent_rank = percent_df(cosmic_rank_dict)
# Aim1.to_csv(cosmic_percent_rank, 'cosmic_percent_rank')






