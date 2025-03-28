import pandas as pd
import pickle
import Aim1
import CheckCanonical as cc
import PhosphositeFeatures as spd

all_mapped_dict = spd.all_mapped_psites
FINAL_sites = cc.FINAL_sites


# function to accept a substrate and a phosphosite and return the sequence and the phosphorylated residue
def get_seq_res(uniprot_id, phosphosite):
    dict_key = str(uniprot_id) + ',' + str(phosphosite)
    if dict_key in all_mapped_dict:
        value = all_mapped_dict.get(dict_key)
        return value


# sequence to find a phosphosite based on a position in the sequence
def get_psite(up_id, position):
    if up_id in FINAL_sites:
        sites_list = FINAL_sites.get(up_id)
        if sites_list is not None:
            for site in sites_list:
                if int(site) - 6 < position < int(site) + 4:
                    return site
    else:
        return None


# method that returns a df of variants so long as they exist at a phosphosite or within +/- 5 AA
def find_vars(df):
    var_df = pd.DataFrame(columns=['UNIPROT_ID', 'SEQUENCE', 'PHOSPHOSITE', 'PS_RESIDUE', 'VAR_POSITION', 'WT_RESIDUE', 'MUT_RESIDUE'])
    #var_df = pd.DataFrame(columns=['UNIPROT_ID', 'SEQUENCE', 'PHOSPHOSITE', 'PS_RESIDUE', 'VAR_POSITION', 'WT_RESIDUE', 'MUT_RESIDUE', 'ENSEMBL_ID', 'TISSUE', 'DISEASE'])
    for item, row in df.iterrows():
        up_id = str(row['UniProt_accession'])
        print(up_id)
        pos = int(row['UniProt_position'])
        wt_res = str(row['residue_WT'])
        mut_res = str(row['residue_variant'])
        # ensembl = row['ENSEMBL_gene_ID']
        # tissue = row['tissue']
        # disease = row['disease_(NCI_thesaurus)']
        psite = get_psite(up_id, pos)
        if psite is not None:
            seq_res = get_seq_res(up_id, psite)
            seq = seq_res[:seq_res.index(',')]
            phos_res = seq_res[seq_res.index(',')+1:]

            data = {
                'UNIPROT_ID': [up_id],
                'SEQUENCE': [seq],
                'PHOSPHOSITE': [psite],
                'PS_RESIDUE': [phos_res],
                'VAR_POSITION': [pos],
                'WT_RESIDUE': [wt_res],
                'MUT_RESIDUE': [mut_res],
                # 'ENSEMBL_ID': [ensembl],
                # 'TISSUE': [tissue],
                # 'DISEASE': [disease]
            }

            df_new_rows = pd.DataFrame(data)
            var_df = pd.concat([var_df, df_new_rows])
    variant_df = var_df.reset_index(drop=True)
    return variant_df


#test_vars = find_vars(cc.header_filtered)
#cosmic = Aim1.cosmic_all_psites
#cosmic_vars_df = find_vars(cosmic)
#Aim1.to_csv(cosmic_vars_df, 'cosmic_variants')

#cosmic2 = Aim1.cosmic_filtered
#cosmic2_vars = find_vars(cosmic2)
#Aim1.to_csv(cosmic2_vars, 'cosmic_variants2')

#cosmic_vars_df = Aim1.cosmic_phos_muts

#cosmic_sample = Aim1.cosmic_samples
#cosmic_sample_vars = find_vars(cosmic_sample)
#Aim1.to_csv(cosmic_sample_vars, 'PM_cosmic_sample')







