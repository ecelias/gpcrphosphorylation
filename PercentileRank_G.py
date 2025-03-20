import pandas as pd
import EvaluateMut_G as em
from scipy.stats import norm
import Aim1_G
import EvaluatePSSM_G as ep
import ScoreDistribution as sd
import time

norm_kinase_list = em.all_kinases
st_norm_kinase_list = ep.st_kinase_list
y_norm_kinase_list = ep.y_kinase_list
sd_st = sd.st_score_dist
sd_y = sd.y_score_dist


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


# function to accept a df with a sequence and PSSM score for each kinase
# returns the percentile for the variant sequence and wt sequence compared to canonical phosphosites
def percent_rank_dict(df):
    percent_dict = {}
    start = time.time()
    for index, row in df.iterrows():
        kinase_dict = {}
        for kinase in kinase_list:
            if kinase in st_kinase_list:
                mean = (sd_st[kinase].astype(float)).mean()
                std = (sd_st[kinase].astype(float)).std()
                score = row[kinase]
                z_score = (score - mean) / std
                percent = norm.cdf(z_score) * 100
                kinase_dict[kinase] = percent
            elif kinase in y_kinase_list:
                mean = (sd_y[kinase].astype(float)).mean()
                std = (sd_y[kinase].astype(float)).std()
                score = row[kinase]
                z_score = (score - mean) / std
                percent = norm.cdf(z_score) * 100
                kinase_dict[kinase] = percent

        key = str(row['ACC_ID']) + ',' + str(row['WT_SEQUENCE']) + ',' + str(row['VARIANT_POSITION']) + ',' + str(
            row['WT_RES']) + ',' + str(row['MUT_RES']) + ',' + str(row['PHOSPHOSITE']) + ',' + str(
            row['PHOSPHOSITE_RESIDUE']) + ',' + str(row['ALLELE_COUNT']) + ',' + str(row['ALLELE_NUM']) + ',' + str(
            row['ALLELE_FREQ']) + ',' + str(row['NUM_HOMOZYG'])
        percent_dict[key] = kinase_dict
        print(row['ACC_ID'])
    end = time.time()
    print(end - start)
    return percent_dict


# transform percent rank dictionary into df
def percent_df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    df = df.reset_index()
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


gnomad_rank_dict = percent_rank_dict(Aim1_G.gnomad_scores)
gnomad_percent_rank = percent_df(gnomad_rank_dict)
Aim1_G.to_csv(gnomad_percent_rank, 'gnomad_percent_rank')
