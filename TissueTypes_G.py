import Aim1
import Aim1_G
import pandas as pd
import pickle

gnomad_phos_muts = Aim1_G.gnomad_phos_muts
gnomad_full = Aim1_G.gnomad_filtered
mapped_tissues = Aim1.tissues_mapped
tpm_unfiltered = Aim1.gtex_tpm
cosmic_tissues = Aim1.cosmic_tissues
tpm = Aim1.tpm_filtered
sty_kinases = Aim1.ser_thr_tyr_kinases


# take a variant df and return the correlating ENSEMBL IDs, tissue type, and disease
def var_tissues(phos_muts_df, variant_df_with_ENSMEBL):
    temp_df = pd.DataFrame(columns=['UNIPROT_ID', 'ENSEMBL_GENE_ID','ALLELE_COUNT', 'ALLELE_NUM', 'ALLELE_FREQ', 'NUM_HOMOZYG'])
    for index, row in phos_muts_df.iterrows():
        current_id = row['UNIPROT_ID']
        count = row['ALLELE_COUNT']
        number = row['ALLELE_NUM']
        freq = row['ALLELE_FREQ']
        homos = row['NUM_HOMOZYG']

        # get ensembl gene id and tissue types
        row = variant_df_with_ENSMEBL.index.get_loc(variant_df_with_ENSMEBL[variant_df_with_ENSMEBL['UniProt_accession'] == current_id].index[0])
        ensembl_id = variant_df_with_ENSMEBL.at[row, 'ENSEMBL_gene_ID']


        data = {
            'UNIPROT_ID': [current_id],
            'ENSEMBL_GENE_ID': [ensembl_id],
            'ALLELE_COUNT': [count],
            'ALLELE_NUM': [number],
            'ALLELE_FREQ': [freq],
            'NUM_HOMOZYG': [homos]
        }

        df_new_rows = pd.DataFrame(data)
        temp_df = pd.concat([temp_df, df_new_rows])
    tissue_df = temp_df.reset_index(drop=True)
    return tissue_df


Aim1.to_csv(var_tissues(gnomad_phos_muts, gnomad_full), 'gnomad_tissues')


# method to get dictionary of cosmic tissue types paired with GTEX tissue types
def tissue_dict():
    tissue_dict = {}
    for index, row in mapped_tissues.iterrows():
        current_tissue = row['COSMIC_tissue']
        if row['GTEx_tissue(s)'] != '-':
            tissue_dict[current_tissue] = row['GTEx_tissue(s)']
    return tissue_dict


# filter gtex tpm
def filter_tpm():
    for index, row in tpm.iterrows():
        ensembl_id = row['Name']
        ensembl_id = ensembl_id[:ensembl_id.index('.')]
        tpm.at[index, 'Name'] = ensembl_id
    return tpm


#mapped_tissue_dict = tissue_dict()
#pickle.dump(mapped_tissue_dict, open("mapped_tissue_dict.p", 'wb'))
cosmic_gtex_mapped_tissue = pickle.load(open("pickle_mapped_tissue_dict.p", "rb"))


def variant_tpm_mapping(tissue_df):
    temp_df2 = pd.DataFrame(
        columns=['UNIPROT_ID', 'ENSEMBL_GENE_ID', 'COSMIC_TISSUE', 'GTEX_TISSUE', 'DISEASE', 'ADIPOSE_TISSUE',
                 'ADRENAL_GLAND', 'BLADDER', 'BLOOD', 'BLOOD_VESSEL', 'BRAIN', 'BREAST', 'CERVIX_UTERI', 'COLON',
                 'ESOPHAGUS', 'FALLOPIAN_TUBE', 'HEART', 'KIDNEY', 'LIVER', 'LUNG', 'MUSCLE', 'NERVE', 'OVARY',
                 'PANCREAS', 'PITUITARY', 'PROSTATE', 'SALIVARY_GLAND', 'SKIN', 'SMALL_INTESTINE', 'SPLEEN', 'STOMACH',
                 'TESTIS', 'THYROID', 'UTERUS', 'VAGINA'])
    for index, row in tissue_df.iterrows():
        current_id = row['UNIPROT_ID']
        ensembl_id = row['ENSEMBL_GENE_ID']
        ensembl_id = ensembl_id[:ensembl_id.index('.')]
        tissue = row['TISSUE']
        disease = row['DISEASE']
        allele_num = row['ALLELE_NUM']
        count = row['ALLELE_COUNT']
        number = row['ALLELE_NUM']
        freq = row['ALLELE_FREQ']
        homos = row['NUM_HOMOZYG']

        gtex_tissue = cosmic_gtex_mapped_tissue.get(tissue)

        col_names = list(tpm['Name'].tolist())

        if ensembl_id in col_names:
            row = tpm.index.get_loc(tpm[tpm['Name'] == ensembl_id].index[0])

            fat = tpm.at[row, 'Adipose Tissue']
            adrenal = tpm.at[row, 'Adrenal Gland']
            bladder = tpm.at[row, 'Bladder']
            blood = tpm.at[row, 'Blood']
            bv = tpm.at[row, 'Blood Vessel']
            brain = tpm.at[row, 'Brain']
            breast = tpm.at[row, 'Breast']
            cervix = tpm.at[row, 'Cervix Uteri']
            colon = tpm.at[row, 'Colon']
            esophagus = tpm.at[row, 'Esophagus']
            fallopian = tpm.at[row, 'Fallopian Tube']
            heart = tpm.at[row, 'Heart']
            kidney = tpm.at[row, 'Kidney']
            liver = tpm.at[row, 'Liver']
            lung = tpm.at[row, 'Lung']
            muscle = tpm.at[row, 'Muscle']
            nerve = tpm.at[row, 'Nerve']
            ovary = tpm.at[row, 'Ovary']
            pancreas = tpm.at[row, 'Pancreas']
            pituitary = tpm.at[row, 'Pituitary']
            prostate = tpm.at[row, 'Prostate']
            saliva = tpm.at[row, 'Salivary Gland']
            skin = tpm.at[row, 'Skin']
            s_intestine = tpm.at[row, 'Small Intestine']
            spleen = tpm.at[row, 'Spleen']
            stomach = tpm.at[row, 'Stomach']
            testis = tpm.at[row, 'Testis']
            thyroid = tpm.at[row, 'Thyroid']
            uterus = tpm.at[row, 'Uterus']
            vagina = tpm.at[row, 'Vagina']

            data = {
                'UNIPROT_ID': [current_id],
                'ENSEMBL_GENE_ID': [ensembl_id],
                'GNOMAD_TISSUE': [tissue],
                'GTEX_TISSUE': [gtex_tissue],
                'DISEASE': [disease],
                'ADIPOSE_TISSUE': [fat],
                'ADRENAL_GLAND': [adrenal],
                'BLADDER': [bladder],
                'BLOOD': [blood],
                'BLOOD_VESSEL': [bv],
                'BRAIN': [brain],
                'BREAST': [breast],
                'CERVIX_UTERI': [cervix],
                'COLON': [colon],
                'ESOPHAGUS': [esophagus],
                'FALLOPIAN_TUBE': [fallopian],
                'HEART': [heart],
                'KIDNEY': [kidney],
                'LIVER': [liver],
                'LUNG': [lung],
                'MUSCLE': [muscle],
                'NERVE': [nerve],
                'OVARY': [ovary],
                'PANCREAS': [pancreas],
                'PITUITARY': [pituitary],
                'PROSTATE': [prostate],
                'SALIVARY_GLAND': [saliva],
                'SKIN': [skin],
                'SMALL_INTESTINE': [s_intestine],
                'SPLEEN': [spleen],
                'STOMACH': [stomach],
                'TESTIS': [testis],
                'THYROID': [thyroid],
                'UTERUS': [uterus],
                'VAGINA': [vagina],
            }

            df_new_rows = pd.DataFrame(data)
            print(df_new_rows)
            temp_df2 = pd.concat([temp_df2, df_new_rows])
    var_tissue_df = temp_df2.reset_index(drop=True)
    return var_tissue_df
