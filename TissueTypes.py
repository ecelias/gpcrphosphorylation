import Aim1
import pandas as pd
import pickle

cosmic_phos_muts = Aim1.cosmic_phos_muts
cosmic_head = Aim1.cosmic_head
cosmic_full = Aim1.cosmic_all_psites
mapped_tissues = Aim1.tissues_mapped
tpm_unfiltered = Aim1.gtex_tpm
cosmic_tissues = Aim1.cosmic_tissues
tpm = Aim1.tpm_filtered
sty_kinases = Aim1.ser_thr_tyr_kinases
all_psites_mapped_ensembl = Aim1.all_psites_mapped_ensembl

sty_tpm = Aim1.sty_kinase_tpm
cosmic_tpm = Aim1.cosmic_tpm
tpm = Aim1.tpm_filtered
all_psites_tpm = Aim1.psite_tpm

# take a variant df and return the correlating ENSEMBL IDs, tissue type, and disease
def var_tissues(phos_muts_df, variant_df_with_ENSMEBL):
    temp_df = pd.DataFrame(columns=['UNIPROT_ID', 'ENSEMBL_GENE_ID', 'TISSUE', 'DISEASE'])
    for index, row in phos_muts_df.iterrows():
        current_id = row['UNIPROT_ID']
        # get ensembl gene id and tissue types
        row = variant_df_with_ENSMEBL.index.get_loc(variant_df_with_ENSMEBL[variant_df_with_ENSMEBL['UniProt_accession'] == current_id].index[0])
        ensembl_id = variant_df_with_ENSMEBL.at[row, 'ENSEMBL_gene_ID']
        tissue = variant_df_with_ENSMEBL.at[row, 'tissue']
        tissue = str(tissue)
        disease = variant_df_with_ENSMEBL.at[row, 'disease_(NCI_thesaurus)']

        data = {
            'UNIPROT_ID': [current_id],
            'ENSEMBL_GENE_ID': [ensembl_id],
            'TISSUE': [tissue],
            'DISEASE': [disease],
        }

        df_new_rows = pd.DataFrame(data)
        temp_df = pd.concat([temp_df, df_new_rows])
    tissue_df = temp_df.reset_index(drop=True)
    return tissue_df


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
                'COSMIC_TISSUE': [tissue],
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


def kinase_tpm_mapping(kinase_tissue_df):
    temp_df3 = pd.DataFrame(
        columns=['KINASE_MATRIX_NAME', 'KINASE_GENE_NAME', 'UNIPROT_ID', 'ENSEMBL_GENE_ID', 'ADIPOSE_TISSUE',
                 'ADRENAL_GLAND', 'BLADDER', 'BLOOD', 'BLOOD_VESSEL', 'BRAIN', 'BREAST', 'CERVIX_UTERI', 'COLON',
                 'ESOPHAGUS', 'FALLOPIAN_TUBE', 'HEART', 'KIDNEY', 'LIVER', 'LUNG', 'MUSCLE', 'NERVE', 'OVARY',
                 'PANCREAS', 'PITUITARY', 'PROSTATE', 'SALIVARY_GLAND', 'SKIN', 'SMALL_INTESTINE', 'SPLEEN', 'STOMACH',
                 'TESTIS', 'THYROID', 'UTERUS', 'VAGINA'])

    for index, row in kinase_tissue_df.iterrows():
        kinase_matrix = row['#Matrix_name']
        print(kinase_matrix)
        kinase_gene = row['Gene_name']
        current_id = row['UniProt']
        ensembl_id = row['Ensembl']

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
                'KINASE_MATRIX_NAME': [kinase_matrix],
                'KINASE_GENE_NAME': [kinase_gene],
                'UNIPROT_ID': [current_id],
                'ENSEMBL_GENE_ID': [ensembl_id],
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
            temp_df3 = pd.concat([temp_df3, df_new_rows])
    kinase_tissue_df = temp_df3.reset_index(drop=True)
    kinase_tissue_df = kinase_tissue_df.set_index('KINASE_MATRIX_NAME')
    return kinase_tissue_df


def psite_tpm_mapping(tissue_df):
    temp_df4 = pd.DataFrame(
        columns=['UNIPROT_ID', 'ENSEMBL_GENE_ID', 'ADIPOSE_TISSUE',
                 'ADRENAL_GLAND', 'BLADDER', 'BLOOD', 'BLOOD_VESSEL', 'BRAIN', 'BREAST', 'CERVIX_UTERI', 'COLON',
                 'ESOPHAGUS', 'FALLOPIAN_TUBE', 'HEART', 'KIDNEY', 'LIVER', 'LUNG', 'MUSCLE', 'NERVE', 'OVARY',
                 'PANCREAS', 'PITUITARY', 'PROSTATE', 'SALIVARY_GLAND', 'SKIN', 'SMALL_INTESTINE', 'SPLEEN', 'STOMACH',
                 'TESTIS', 'THYROID', 'UTERUS', 'VAGINA'])
    for index, row in tissue_df.iterrows():
        current_id = row['UNIPROT_ID']
        ensembl_id = row['ENSEMBL_ID']

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
            temp_df4 = pd.concat([temp_df4, df_new_rows])
    psite_tissue_df = temp_df4.reset_index(drop=True)
    psite_tissue_df = psite_tissue_df.set_index('UNIPROT_ID')
    return psite_tissue_df


#Aim1.to_csv(kinase_tpm_mapping(sty_kinases), 'sty_kinase_tpm_mapped')
#Aim1.to_csv(psite_tpm_mapping(all_psites_mapped_ensembl), 'psite_tpm_mapped')

def proteins_of_interest(psite_tpms, kinase_tpms):
    df = pd.DataFrame(
        columns=['UNIPROT_ID', 'ENSEMBL_ID', 'ADIPOSE_TISSUE', 'ADRENAL_GLAND', 'BLADDER', 'BLOOD', 'BLOOD_VESSEL',
                 'BRAIN', 'BREAST', 'CERVIX_UTERI', 'COLON', 'ESOPHAGUS', 'FALLOPIAN_TUBE', 'HEART', 'KIDNEY', 'LIVER',
                 'LUNG', 'MUSCLE', 'NERVE', 'OVARY', 'PANCREAS', 'PITUITARY', 'PROSTATE', 'SALIVARY_GLAND', 'SKIN',
                 'SMALL_INTESTINE', 'SPLEEN', 'STOMACH', 'TESTIS', 'THYROID', 'UTERUS', 'VAGINA'])
    for index, row in psite_tpms.iterrows():
        up_id = row['UNIPROT_ID']
        ensembl_id = row['ENSEMBL_GENE_ID']
        print(ensembl_id)
        r = tpm.index.get_loc(tpm[tpm['Name'] == ensembl_id].index[0])

        fat = tpm.at[r, 'Adipose Tissue']
        adrenal = tpm.at[r, 'Adrenal Gland']
        bladder = tpm.at[r, 'Bladder']
        blood = tpm.at[r, 'Blood']
        bv = tpm.at[r, 'Blood Vessel']
        brain = tpm.at[r, 'Brain']
        breast = tpm.at[r, 'Breast']
        cervix = tpm.at[r, 'Cervix Uteri']
        colon = tpm.at[r, 'Colon']
        esophagus = tpm.at[r, 'Esophagus']
        fallopian = tpm.at[r, 'Fallopian Tube']
        heart = tpm.at[r, 'Heart']
        kidney = tpm.at[r, 'Kidney']
        liver = tpm.at[r, 'Liver']
        lung = tpm.at[r, 'Lung']
        muscle = tpm.at[r, 'Muscle']
        nerve = tpm.at[r, 'Nerve']
        ovary = tpm.at[r, 'Ovary']
        pancreas = tpm.at[r, 'Pancreas']
        pituitary = tpm.at[r, 'Pituitary']
        prostate = tpm.at[r, 'Prostate']
        saliva = tpm.at[r, 'Salivary Gland']
        skin = tpm.at[r, 'Skin']
        s_intestine = tpm.at[r, 'Small Intestine']
        spleen = tpm.at[r, 'Spleen']
        stomach = tpm.at[r, 'Stomach']
        testis = tpm.at[r, 'Testis']
        thyroid = tpm.at[r, 'Thyroid']
        uterus = tpm.at[r, 'Uterus']
        vagina = tpm.at[r, 'Vagina']

        data = {
            'UNIPROT_ID': [up_id],
            'ENSEMBL_ID': [ensembl_id],
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
        df = pd.concat([df, df_new_rows])

    for index, row in kinase_tpms.iterrows():
        up_id = row['UNIPROT_ID']
        ensembl_id_kin = row['ENSEMBL_GENE_ID']
        print(ensembl_id_kin)

        r = tpm.index.get_loc(tpm[tpm['Name'] == ensembl_id_kin].index[0])

        fat = tpm.at[r, 'Adipose Tissue']
        adrenal = tpm.at[r, 'Adrenal Gland']
        bladder = tpm.at[r, 'Bladder']
        blood = tpm.at[r, 'Blood']
        bv = tpm.at[r, 'Blood Vessel']
        brain = tpm.at[r, 'Brain']
        breast = tpm.at[r, 'Breast']
        cervix = tpm.at[r, 'Cervix Uteri']
        colon = tpm.at[r, 'Colon']
        esophagus = tpm.at[r, 'Esophagus']
        fallopian = tpm.at[r, 'Fallopian Tube']
        heart = tpm.at[r, 'Heart']
        kidney = tpm.at[r, 'Kidney']
        liver = tpm.at[r, 'Liver']
        lung = tpm.at[r, 'Lung']
        muscle = tpm.at[r, 'Muscle']
        nerve = tpm.at[r, 'Nerve']
        ovary = tpm.at[r, 'Ovary']
        pancreas = tpm.at[r, 'Pancreas']
        pituitary = tpm.at[r, 'Pituitary']
        prostate = tpm.at[r, 'Prostate']
        saliva = tpm.at[r, 'Salivary Gland']
        skin = tpm.at[r, 'Skin']
        s_intestine = tpm.at[r, 'Small Intestine']
        spleen = tpm.at[r, 'Spleen']
        stomach = tpm.at[r, 'Stomach']
        testis = tpm.at[r, 'Testis']
        thyroid = tpm.at[r, 'Thyroid']
        uterus = tpm.at[r, 'Uterus']
        vagina = tpm.at[r, 'Vagina']

        data = {
            'UNIPROT_ID': [up_id],
            'ENSEMBL_ID': [ensembl_id_kin],
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
        df = pd.concat([df, df_new_rows])
    poi_df = df.reset_index(drop=True)
    return poi_df


#Aim1.to_csv(proteins_of_interest(Aim1.psite_tpm, sty_tpm), 'proteins_of_interest')