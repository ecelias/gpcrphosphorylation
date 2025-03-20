import pandas as pd
import Aim1
import EvaluateMut as em
import EvaluatePSSM as ep

pr_1 = Aim1.readCSV("/home/eelias13/KinaseModulation/PR_updated_files/cosmic_percent_rank_with_updated_scores_0_to_50.csv")
pr_2 = Aim1.readCSV("/home/eelias13/KinaseModulation/PR_updated_files/cosmic_percent_rank_with_updated_scores_50_to_100.csv")
pr_3 = Aim1.readCSV("/home/eelias13/KinaseModulation/PR_updated_files/cosmic_percent_rank_with_updated_scores_100_to_150.csv")
pr_4 = Aim1.readCSV("/home/eelias13/KinaseModulation/PR_updated_files/cosmic_percent_rank_with_updated_scores_150_to_200.csv")

updated_pr = pd.concat([pr_1, pr_2, pr_3, pr_4])
norm_kinase_list = em.all_kinases
st_kinase_list = ep.st_kinase_list
y_kinase_list = ep.y_kinase_list

def mutant_features(index_key):
    features = index_key.split('_')
    return features

def percent_change(df, threshold):
    change_dict = {}
    phos_mim_or_null = pd.DataFrame()
    for index, row in df.iterrows():
        kinase_dict = {}

        # mut_features = mutant_features(row['Unnamed: 0'])
        # up_id = mut_features[0]
        # seq = mut_features[1]
        # var_info = mut_features[4]
        # var_pos = var_info[1:-1]
        # wt_residue = var_info[:1]
        # mut_residue = var_info[-1:]

        # psite = mut_features[2]
        # p_residue = mut_features[3]

        up_id = row['ACC_ID']
        seq =  row['WT_SEQUENCE']
        var_pos = row['VARIANT_POSITION']
        wt_residue = row['WT_RES']
        mut_residue = row['MUT_RES']

        psite = row['PHOSPHOSITE']
        p_residue =row['PHOSPHOSITE_RESIDUE']

        if var_pos == psite and mut_residue not in {'Y', 'S', 'T'}:
            data = {
                'UNIPROT_ID': [up_id],
                'SEQUENCE': [seq],
                'PHOSPHOSITE': [psite],
                'PS_RESIDUE': [p_residue],
                'VAR_POSITION': [var_pos],
                'WT_RESIDUE': [wt_residue],
                'MUT_RESIDUE': [mut_residue],
            }
            df_new_rows = pd.DataFrame(data)
            phos_mim_or_null = pd.concat([phos_mim_or_null, df_new_rows])
        else:
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]

                wt_score = wt_score if pd.notnull(wt_score) else 0.0
                mut_score = mut_score if pd.notnull(mut_score) else 0.0

                # only count if either mut or wt above threshold, otherwise set to 0
                if wt_score >= float(threshold) and mut_score <= float(threshold):
                    change = wt_score - threshold
                    # change between threshold and number
                    kinase_dict[kinase] = change
                elif mut_score >= float(threshold) and wt_score <= float(threshold):
                    change = mut_score - threshold
                    # change between threshold and number
                    kinase_dict[kinase] = change
                elif wt_score > float(threshold) and mut_score > float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change
                elif wt_score < float(threshold) and mut_score < float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change

            k = str(up_id) + '_' + seq + '_' + str(psite) + '_' + p_residue + '_' + wt_residue + str(var_pos) + mut_residue
            change_dict[k] = kinase_dict
    phos_mim_or_null = phos_mim_or_null.reset_index(drop=True)
    print(phos_mim_or_null)
    return change_dict

def percent_change_phosm_null(df, threshold):
    change_dict = {}
    phos_mim_or_null = pd.DataFrame()
    for index, row in df.iterrows():
        kinase_dict = {}

        # mut_features = mutant_features(row['Unnamed: 0'])
        # up_id = mut_features[0]
        # seq = mut_features[1]
        # var_info = mut_features[4]
        # var_pos = var_info[1:-1]
        # wt_residue = var_info[:1]
        # mut_residue = var_info[-1:]

        # psite = mut_features[2]
        # p_residue = mut_features[3]

        up_id = row['ACC_ID']
        seq =  row['WT_SEQUENCE']
        var_pos = row['VARIANT_POSITION']
        wt_residue = row['WT_RES']
        mut_residue = row['MUT_RES']

        psite = row['PHOSPHOSITE']
        p_residue =row['PHOSPHOSITE_RESIDUE']

        if var_pos == psite and mut_residue not in {'Y', 'S', 'T'}:
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]

                wt_score = wt_score if pd.notnull(wt_score) else 0.0
                mut_score = mut_score if pd.notnull(mut_score) else 0.0


                if p_residue in {'S', 'T'}:
                    # only count if either mut or wt above threshold, otherwise set to 0
                    if wt_score >= float(threshold):
                        if kinase in st_kinase_list:
                            change = wt_score - threshold
                            # change between threshold and number
                            kinase_dict[kinase] = change
                        else: 
                            change = 0.0
                            # change between threshold and number
                            kinase_dict[kinase] = change
                    else: 
                        change = 0.0
                        kinase_dict[kinase] = change

                elif p_residue in {'Y'}:
                    # only count if either mut or wt above threshold, otherwise set to 0
                    if wt_score >= float(threshold):
                        if kinase in y_kinase_list:
                            change = wt_score - threshold
                            # change between threshold and number
                            kinase_dict[kinase] = change
                        else: 
                            change = 0.0
                            # change between threshold and number
                            kinase_dict[kinase] = change
                    else: 
                        change = 0.0
                        kinase_dict[kinase] = change

            k = str(up_id) + '_' + seq + '_' + str(psite) + '_' + p_residue + '_' + wt_residue + str(var_pos) + mut_residue
            change_dict[k] = kinase_dict
    return change_dict

def percent_change_flanking_ST(df, threshold):
    change_dict = {}
    for index, row in df.iterrows():
        kinase_dict = {}

        # mut_features = mutant_features(row['Unnamed: 0'])
        # var_info = mut_features[4]
        # var_pos = var_info[1:-1]
        # wt_residue = var_info[:1]
        # mut_residue = var_info[-1:]

        # psite = mut_features[2]
        # p_residue = mut_features[3]
        # up_id = mut_features[0]
        # seq = mut_features[1]
        up_id = row['ACC_ID']
        seq =  row['WT_SEQUENCE']
        var_pos = row['VARIANT_POSITION']
        wt_residue = row['WT_RES']
        mut_residue = row['MUT_RES']
        print(up_id)

        psite = row['PHOSPHOSITE']
        p_residue =row['PHOSPHOSITE_RESIDUE']

        if var_pos != psite and p_residue in {'S', 'T'}:
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]

                wt_score = wt_score if pd.notnull(wt_score) else 0.0
                mut_score = mut_score if pd.notnull(mut_score) else 0.0

                # only count if either mut or wt above threshold, otherwise set to 0
                if wt_score >= float(threshold) and mut_score <= float(threshold):
                    if kinase in st_kinase_list:
                        change = threshold - wt_score
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else: 
                        change = 0.0
                        # change between threshold and number
                        kinase_dict[kinase] = change
                elif mut_score >= float(threshold) and wt_score <= float(threshold):
                    if kinase in st_kinase_list:
                        change = mut_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else: 
                        change = 0.0
                        # change between threshold and number
                        kinase_dict[kinase] = change 
                elif wt_score > float(threshold) and mut_score > float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change
                elif wt_score < float(threshold) and mut_score < float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change

            k = str(up_id) + '_' + seq + '_' + str(psite) + '_' + p_residue + '_' + wt_residue + str(var_pos) + mut_residue
            change_dict[k] = kinase_dict
    return change_dict


def percent_change_flanking_Y(df, threshold):
    change_dict = {}
    for index, row in df.iterrows():
        kinase_dict = {}

        # mut_features = mutant_features(row['Unnamed: 0'])
        # var_info = mut_features[4]
        # var_pos = var_info[1:-1]
        # wt_residue = var_info[:1]
        # mut_residue = var_info[-1:]

        # psite = mut_features[2]
        # p_residue = mut_features[3]
        # up_id = mut_features[0]
        # seq = mut_features[1]
        up_id = row['ACC_ID']
        seq =  row['WT_SEQUENCE']
        var_pos = row['VARIANT_POSITION']
        wt_residue = row['WT_RES']
        mut_residue = row['MUT_RES']

        psite = row['PHOSPHOSITE']
        p_residue =row['PHOSPHOSITE_RESIDUE']

        if var_pos != psite and p_residue == 'Y':
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]

                wt_score = wt_score if pd.notnull(wt_score) else 0.0
                mut_score = mut_score if pd.notnull(mut_score) else 0.0

                # only count if either mut or wt above threshold, otherwise set to 0
                if wt_score >= float(threshold) and mut_score <= float(threshold):
                    if kinase in y_kinase_list:
                        change = wt_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else: 
                        change = 0.0
                        # change between threshold and number
                        kinase_dict[kinase] = change
                elif mut_score >= float(threshold) and wt_score <= float(threshold):
                    if kinase in y_kinase_list:
                        change = mut_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else: 
                        change = 0.0
                        # change between threshold and number
                        kinase_dict[kinase] = change 
                elif wt_score > float(threshold) and mut_score > float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change
                elif wt_score < float(threshold) and mut_score < float(threshold):
                    change = 0.0
                    kinase_dict[kinase] = change

            k = str(up_id) + '_' + seq + '_' + str(psite) + '_' + p_residue + '_' + wt_residue + str(var_pos) + mut_residue
            change_dict[k] = kinase_dict
    return change_dict

def percent_change_st_to_y(df, threshold):
    change_dict = {}
    for index, row in df.iterrows():
        kinase_dict = {}

        # mut_features = mutant_features(row['Unnamed: 0'])
        # var_info = mut_features[4]
        # var_pos = var_info[1:-1]
        # wt_residue = var_info[:1]
        # mut_residue = var_info[-1:]

        # psite = mut_features[2]
        # p_residue = mut_features[3]
        # up_id = mut_features[0]
        # seq = mut_features[1]
        up_id = row['ACC_ID']
        seq =  row['WT_SEQUENCE']
        var_pos = row['VARIANT_POSITION']
        wt_residue = row['WT_RES']
        mut_residue = row['MUT_RES']

        psite = row['PHOSPHOSITE']
        p_residue =row['PHOSPHOSITE_RESIDUE']

        if var_pos == psite and wt_residue in {'S', 'T'} and mut_residue == 'Y':
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]

                wt_score = wt_score if pd.notnull(wt_score) else 0.0
                mut_score = mut_score if pd.notnull(mut_score) else 0.0

                # only count if either mut or wt above threshold, otherwise set to 0
                if wt_score >= float(threshold):
                    if kinase in st_kinase_list:
                        change = wt_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else:
                        change = 0.0
                        kinase_dict[kinase] = change
                elif mut_score >= float(threshold):
                    if kinase in y_kinase_list:
                        change = mut_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else:
                        change = 0.0
                        kinase_dict[kinase] = change
                else:
                    change = 0.0
                    kinase_dict[kinase] = change

            k = str(up_id) + '_' + seq + '_' + str(psite) + '_' + p_residue + '_' + wt_residue + str(var_pos) + mut_residue
            change_dict[k] = kinase_dict
    return change_dict


def percent_change_y_to_st(df, threshold):
    change_dict = {}
    for index, row in df.iterrows():
        kinase_dict = {}

        #mut_features = mutant_features(row['Unnamed: 0'])
        # var_info = mut_features[4]
        # var_pos = var_info[1:-1]
        # wt_residue = var_info[:1]
        # mut_residue = var_info[-1:]

        # psite = mut_features[2]
        # p_residue = mut_features[3]
        # up_id = mut_features[0]
        # seq = mut_features[1]

        up_id = row['ACC_ID']
        seq =  row['WT_SEQUENCE']
        var_pos = row['VARIANT_POSITION']
        wt_residue = row['WT_RES']
        mut_residue = row['MUT_RES']

        psite = row['PHOSPHOSITE']
        p_residue =row['PHOSPHOSITE_RESIDUE']

        if var_pos == psite and mut_residue in {'S', 'T'} and wt_residue == 'Y':
            for kinase in norm_kinase_list:
                mut_kinase = kinase + '_m'
                wt_score = row[kinase]
                mut_score = row[mut_kinase]

                wt_score = wt_score if pd.notnull(wt_score) else 0.0
                mut_score = mut_score if pd.notnull(mut_score) else 0.0

                # only count if either mut or wt above threshold, otherwise set to 0
                if wt_score >= float(threshold):
                    if kinase in y_kinase_list:
                        change = wt_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else:
                        change = 0.0
                        kinase_dict[kinase] = change
                elif mut_score >= float(threshold):
                    if kinase in st_kinase_list:
                        change = mut_score - threshold
                        # change between threshold and number
                        kinase_dict[kinase] = change
                    else:
                        change = 0.0
                        kinase_dict[kinase] = change
                else:
                    change = 0.0
                    kinase_dict[kinase] = change

            k = str(up_id) + '_' + seq + '_' + str(psite) + '_' + p_residue + '_' + wt_residue + str(var_pos) + mut_residue
            change_dict[k] = kinase_dict
    return change_dict


def change_df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    pd.options.display.float_format = '{:.3f}'.format
    #df2 = df.reset_index(drop=True)
    #print(df2)
    #return df2
    return df


# pr_sample = Aim1.readCSV('cosmic_sample_percent_rank.csv')
pr = pd.read_csv('cosmic_percent_rank.csv', index_col=0)

# dv_df_phosm_null = change_df(percent_change_phosm_null(pr, 90))
# dv_df_flank_y = change_df(percent_change_flanking_Y(pr, 90))
# dv_df_flank_st = change_df(percent_change_flanking_ST(pr, 90))
# dv_df_st_to_y = change_df(percent_change_st_to_y(pr, 90))
# dv_df_y_to_st = change_df(percent_change_y_to_st(pr, 90))

# Aim1.to_csv(dv_df_phosm_null, 'dv_cosmic_phosm_null_threshold_90')
# print("done")
# Aim1.to_csv(dv_df_flank_y, 'dv_cosmic_flank_y_threshold_90')
# print("done")
# Aim1.to_csv(dv_df_flank_st, 'dv_cosmic_flank_st_threshold_90')
# print("done")
# Aim1.to_csv(dv_df_st_to_y, 'dv_cosmic_st_to_y_threshold_90')
# print("done")
# Aim1.to_csv(dv_df_y_to_st, 'dv_cosmic_y_to_st_threshold_90')
# print("done")
