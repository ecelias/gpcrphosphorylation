import pandas as pd
import Aim1
import PhosphositeFeatures as spd
import pickle
import time

cosmic_head = Aim1.cosmic_head
cosmic = Aim1.cosmic_full


# create a dictionary for phosphosites where key is a kinase and value is a list of all phosphosites
# dictionary entry: {'Q9ULL0': [1022, 1086, 1441, 1522, 472, 543, 816, 818, 915, 965]}
def phosphosite_dict(df):
    phosphosite_dict = {}
    for item, row in df.iterrows():
        site_id = row['UniProt_accession_mapped']
        phosphosite = row['position_mapped']
        if site_id not in phosphosite_dict:
            phosphosite_dict[site_id] = [phosphosite]
        else:
            phosphosite_dict[site_id].append(phosphosite)
    return phosphosite_dict


# filter the data frame so you only include sequences for the canonical isoform in SwissProt
def filter_vars(df, dict1):
    filtered = pd.DataFrame(columns=df.columns)
    i = 0
    for index, row in df.iterrows():
        current_id = row['UniProt_accession']
        var_position = row['UniProt_position']
        if Aim1.getSeq(current_id) is not None:
            if current_id in dict1.keys():
                sites_list = dict1.get(current_id)
                for site in sites_list:
                    if int(site) - 6 < var_position < int(site) + 4:
                        filtered.loc[i] = df.loc[index].copy()
        i += 1
    filtered = filtered.reset_index(drop=True)
    return filtered


#final_phos_sites = phosphosite_dict(Aim1.final_phosphosites)
#pickle.dump(final_phos_sites, open("pickle_final_phos_sites.p", 'wb'))
FINAL_sites = pickle.load(open("pickle_final_phos_sites.p", "rb"))


#header_filtered = filter_vars(cosmic_head, FINAL_sites)
#filtered_vars_all_psites = filter_vars(cosmic, FINAL_sites)
#Aim1.to_csv(filtered_vars_all_psites, 'cosmic_filtered_all_psites')