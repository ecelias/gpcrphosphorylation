import pandas as pd
import Aim1
import pickle
import functools
import time
import DifferenceVector as dv


tpm = Aim1.tpm_filtered
sty_tpm = Aim1.sty_kinase_tpm
all_psites_tpm = Aim1.psite_tpm

cosmic_pr = Aim1.cosmic_percent_rank
poi = Aim1.poi

gtex_tissue_names = pickle.load(open("pickle_gtex_tissue_names.p", "rb"))
matrix_to_up_ids = pickle.load(open("pickle_matrix_to_up.p", "rb"))
cosmic_up_ids = pickle.load(open("pickle_cosmic_up_ids.p", "rb"))


def poi_dict(df):
    poi_dict = {}
    for index, row in df.iterrows():
        current_id = row['UNIPROT_ID']
        for tissue in gtex_tissue_names:
            current_tpm = float(row[tissue])
            if current_tpm is None:
                current_tpm = 0.0
            if current_id in poi_dict:
                poi_dict[current_id].append(current_tpm)
            else:
                poi_dict[current_id] = [current_tpm]
    return poi_dict


#temp = poi_dict(poi)
#pickle.dump(temp, open("pickle_poi_dict.p", 'wb'))
poi_dic = pickle.load(open("pickle_poi_dict.p", "rb"))


def find_relevance(percent_rank_df, dv_threshold, tpm_threshold):
    start = time.time()

    dv_dict = dv.percent_change(percent_rank_df, dv_threshold)
    dv_df = dv.change_df(dv_dict)
    dv_cleaned = dv_df.drop(columns=[col for col in dv_df if (dv_df[col] == 0).all()])

    for index, row in dv_cleaned.iterrows():
        substrate_id = index.split(',')[0]
        print(substrate_id)
        for kinase in matrix_to_up_ids.keys():
            print(kinase)
            kinase_id = matrix_to_up_ids.get(kinase)

            substrate_tpm = poi_dic.get(substrate_id)
            kinase_tpm = poi_dic.get(kinase_id)
            if kinase_tpm is None or substrate_tpm is None or len(kinase_tpm) != len(substrate_tpm):
                dv_df.loc[index, kinase] = 0
            else:
                for i in range(len(substrate_tpm)):
                    difference_vector = dv_cleaned.loc[index, kinase]
                    if substrate_tpm[i] >= float(tpm_threshold) and kinase_tpm[i] >= float(tpm_threshold):
                        dv_cleaned.loc[index, kinase] = difference_vector
                        break
                    else:
                        dv_df.loc[index, kinase] = 0

    end = time.time()
    elapsed = end - start
    print(elapsed)
    dv_cleaned2 = dv_cleaned.reset_index(drop=True)
    return dv_cleaned2


# function to call the find_relevance function
def get_relevance(percent_rank_df):
    difference_threshold = input('Set threshold for difference vector: ')
    print(f'Difference vector threshold set to {difference_threshold}')

    tpm_threshold = input('Set threshold for tpm: ')
    print(f'TPM threshold set to {tpm_threshold}')

    csv_name = input('What do you want to name the CSV file')

    return Aim1.to_csv(find_relevance(percent_rank_df, difference_threshold, tpm_threshold), csv_name)


Aim1.to_csv(find_relevance(cosmic_pr, 90, 20), 'tissue_relevance_20')
