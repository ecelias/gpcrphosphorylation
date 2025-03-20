import pandas as pd
import Aim1_G
import CheckCanonical_G as ccg
import SubstratePhosphositeDictionary as spd


all_mapped_dict = spd.all_mapped_psites
FINAL_sites = ccg.FINAL_sites


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


def find_vars(df):
    var_df = pd.DataFrame(columns=['UNIPROT_ID', 'SEQUENCE', 'PHOSPHOSITE', 'PS_RESIDUE', 'VAR_POSITION', 'WT_RESIDUE', 'MUT_RESIDUE', 'ALLELE_COUNT', 'ALLELE_NUM', 'ALLELE_FREQ', 'NUM_HOMOZYG'])
    for item, row in df.iterrows():
        up_id = str(row['UniProt_accession'])
        print(up_id)
        pos = int(row['UniProt_position'])
        wt_res = str(row['residue_WT'])
        mut_res = str(row['residue_variant'])
        allele_count = row['allele_count']
        allele_num = row['allele_number']
        allele_freq = row['allele_frequency']
        num_homo = row['number_of_homozygotes']
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
                'ALLELE_COUNT': [allele_count],
                'ALLELE_NUM': [allele_num],
                'ALLELE_FREQ': [allele_freq],
                'NUM_HOMOZYG': [num_homo]
            }
            df_new_rows = pd.DataFrame(data)
            var_df = pd.concat([var_df, df_new_rows])
    variant_df = var_df.reset_index(drop=True)
    return variant_df


#test_vars = find_vars(Aim2_G.filtered_gnomad_sample)

#gnomad_vars_df = find_vars(Aim1_G.gnomad_filtered)
#Aim1_G.to_csv(gnomad_vars_df, 'gnomad_variants')







