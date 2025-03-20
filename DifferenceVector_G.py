import PercentileRank as pr
import pandas as pd
import Aim1
import EvaluateMut as em

kinase_list = pr.kinase_list
test_cosmic_pr = Aim1.readCSV('C:/Users/eelias13/PycharmProjects/KinaseNetworkModulation/pr.test.csv')
norm_kinase_list = em.all_kinases


def percent_change(df, threshold):
    change_dict = {}
    phos_mim_or_null = pd.DataFrame()
    for index, row in df.iterrows():
        kinase_dict = {}
        mut_residue = row['MUT_RES']
        var_pos = row['VARIANT_POSITION']
        psite = row['PHOSPHOSITE']

        if var_pos == psite and mut_residue not in {'Y', 'S', 'T'}:
            data = {
                'UNIPROT_ID': [row['ACC_ID']],
                'SEQUENCE': [row['WT_SEQUENCE']],
                'PHOSPHOSITE': [row['PHOSPHOSITE']],
                'PS_RESIDUE': [row['PHOSPHOSITE_RESIDUE']],
                'VAR_POSITION': [row['VARIANT_POSITION']],
                'WT_RESIDUE': [row['WT_RES']],
                'MUT_RESIDUE': [row['MUT_RES']],
                'ALLELE_COUNT': [row['ALLELE_COUNT']],
                'ALLELE_NUM': [row['ALLELE_NUM']],
                'ALLELE_FREQ': [row['ALLELE_FREQ']],
                'NUM_HOMOZYG': [row['NUM_HOMOZYG']]
                ,
            }
            df_new_rows = pd.DataFrame(data)
            phos_mim_or_null = pd.concat([phos_mim_or_null, df_new_rows])
        else:
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]
                if wt_score >= float(threshold) or mut_score >= float(threshold):
                    change = wt_score - mut_score
                    kinase_dict[kinase] = change
                elif wt_score < float(threshold) and mut_score < float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change
            k = row['WT_SEQUENCE']
            change_dict[k] = kinase_dict
    phos_mim_or_null = phos_mim_or_null.reset_index(drop=True)
    print(phos_mim_or_null)
    return change_dict


def change_df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    pd.options.display.float_format = '{:.3f}'.format
    return df

