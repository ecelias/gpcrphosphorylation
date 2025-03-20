import pandas as pd

import Aim1
import Aim1_G

gnomad = Aim1_G.gnomad
gnomad_filtered = Aim1_G.gnomad_filtered
filtered_vars_g = Aim1_G.gnomad_phos_muts
proteome = Aim1.proteome


# number of mutations in full human proteome
def total_mutations(mutation_df):
    n_psite_mutations = 0
    for index, row in mutation_df.iterrows():
        up_id = row['UniProt_accession']
        print(up_id)
        current_af = row['allele_frequency']
        if Aim1.getSeq(up_id) is not None:
            n_psite_mutations += current_af
    return n_psite_mutations


# number of mutations in the full human proteome that map to a phosphosite
def total_psite_mutations(filtered_mutation_df):
    n_psite_mutations = 0
    for index, row in filtered_mutation_df.iterrows():
        up_id = row['UniProt_accession']
        current_af = row['allele_frequency']
        current_af = row['ALLELE_FREQ']
        if Aim1.getSeq(up_id) is not None:
            n_psite_mutations += current_af
    return n_psite_mutations


# number of null mutations
def null_mutations(filtered_vars_df):
    n_null = 0
    for index, row in filtered_vars_df.iterrows():
        psite = row['PHOSPHOSITE']
        var_pos = row['VAR_POSITION']
        mut_res = row['MUT_RESIDUE']
        current_af = row['ALLELE_FREQ']
        if psite == var_pos and mut_res not in {'Y', 'S', 'D', 'E', 'T'}:
            n_null += current_af
    return n_null


# number of phosphomimetic mutations
def phos_m_mutations(filtered_vars_df):
    n_phos_m = 0
    for index, row in filtered_vars_df.iterrows():
        psite = row['PHOSPHOSITE']
        var_pos = row['VAR_POSITION']
        mut_res = row['MUT_RESIDUE']
        current_af = row['ALLELE_FREQ']
        if psite == var_pos and mut_res not in {'Y', 'S', 'T'} and mut_res in {'D', 'E'}:
            n_phos_m += current_af
    return n_phos_m


# number of mutations that swap S/T -> Y
def to_tyr(filtered_vars_df):
    n_to_tyr = 0
    for index, row in filtered_vars_df.iterrows():
        psite = row['PHOSPHOSITE']
        var_pos = row['VAR_POSITION']
        wt_res = row['WT_RESIDUE']
        mut_res = row['MUT_RESIDUE']
        current_af = row['ALLELE_FREQ']
        if psite == var_pos and wt_res in {'S', 'T'} and mut_res == 'Y':
            n_to_tyr += current_af
    return n_to_tyr


# number of mutations that swap S -> T or T -> S
def ser_to_thr(filtered_vars_df):
    n_to_ser_thr = 0
    for index, row in filtered_vars_df.iterrows():
        psite = row['PHOSPHOSITE']
        var_pos = row['VAR_POSITION']
        wt_res = row['WT_RESIDUE']
        mut_res = row['MUT_RESIDUE']
        current_af = row['ALLELE_FREQ']
        if psite == var_pos and mut_res == 'S' and wt_res == 'T':
            n_to_ser_thr += current_af
        if psite == var_pos and mut_res == 'T' and wt_res == 'S':
            n_to_ser_thr += current_af
    return n_to_ser_thr


# number of mutations that swap Y -> S/T
def to_ser_thr(filtered_vars_df):
    n_to_ser_thr = 0
    for index, row in filtered_vars_df.iterrows():
        psite = row['PHOSPHOSITE']
        var_pos = row['VAR_POSITION']
        wt_res = row['WT_RESIDUE']
        mut_res = row['MUT_RESIDUE']
        current_af = row['ALLELE_FREQ']
        if psite == var_pos and mut_res in {'S', 'T'} and wt_res == 'Y':
            n_to_ser_thr += current_af
    return n_to_ser_thr


# number of mutations in the flanking region
def flanking_mutations(filtered_vars_df):
    n_flank_mut = 0
    for index, row in filtered_vars_df.iterrows():
        psite = row['PHOSPHOSITE']
        var_pos = row['VAR_POSITION']
        current_af = row['ALLELE_FREQ']
        if psite != var_pos:
            n_flank_mut += current_af
    return n_flank_mut


def summary_stats(mutation_df, filtered_mutation_df, filtered_vars_df):
    sum_df = pd.DataFrame(columns=['TOTAL_MUTS', 'TOTAL_PSITE_MUTS', 'NULL_MUTS', 'MIMIC_MUTS', 'ST_to_Y', 'Y_to_ST', 'FLANK_MUTS'])
    data = {
        'TOTAL_MUTS' : [total_mutations(mutation_df)],
        'TOTAL_PSITE_MUTS': [total_psite_mutations(filtered_mutation_df)],
        'NULL_MUTS': [null_mutations(filtered_vars_df)],
        'MIMIC_MUTS': [phos_m_mutations(filtered_vars_df)],
        'ST_to_Y': [to_tyr(filtered_vars_df)],
        'Y_to_ST': [to_ser_thr(filtered_vars_df)],
        'FLANK_MUTS': [flanking_mutations(filtered_vars_df)],
    }
    df_new_rows = pd.DataFrame(data)
    sum_df = pd.concat([sum_df, df_new_rows])
    summary_df = sum_df.reset_index(drop=True)
    return summary_df


#Aim1.to_csv(summary_stats(gnomad, gnomad_filtered, filtered_vars_g), 'gnomad_summary_stats')

print(ser_to_thr(filtered_vars_g))
