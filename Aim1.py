import re
import pandas as pd
import pickle


# FASTA to data frame
# from fastaframes import to_df

#fasta_df = to_df(data="C:/Users/eelias13/PycharmProjects/KinaseNetworkModulation/human_proteome.fasta")


# CSV reader (CSV to data frame)
def readCSV(csvFileName):
    try:
        df = pd.read_csv(csvFileName)
        return df
    except FileNotFoundError:
        print(f"The file {csvFileName} was not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []


# TSV reader (TSV to data frame)
def readTSV(tsvFileName):
    try:
        df = pd.read_csv(tsvFileName, sep='\t')
        return df
    except FileNotFoundError:
        print(f"The file {tsvFileName} was not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []


# csv files
proteomeCSV = "/home/eelias13/KinaseModulation/ProteomeFiles/proteome.csv"

pssmCSV = "/home/eelias13/KinaseModulation/pssmFiles/pssm.csv"
tyrPSSMCSV = '/home/eelias13/KinaseModulation/pssmFiles/tyrosine_pssm.csv'

final_phosphositesCSV = "C:/Users/eelias13/Desktop/KinaseModulation/Mapped_Phosphosites/final_mapped_phosphosites.csv"
all_psites_mapped_ensemblCSV = "C:/Users/eelias13/Desktop/KinaseModulation/Mapped_Phosphosites/all_phosphosites_mapped_to_ensembl_ids.csv"

#cosmicTSV = 'C:/Users/eelias13/PycharmProjects/KinaseModulation/CosmicFiles/COSMIC_UniProt_mapped_GenomeScreen_missense.tsv'
cosmicTSV = ''
cosmic_headCSV = '/home/eelias13/KinaseModulation/CosmicFiles/COSMIC_UniProt_mapped_GenomeScreen_missense-sample.csv'
cosmic_phos_mutsCSV = '/home/eelias13/KinaseModulation/CosmicFiles/cosmic_variants.csv'
cosmic_scoresCSV = 'CosmicFiles/cosmic_kinase_scores.csv'
cosmic_all_psitesCSV = 'C:/Users/eelias13/PycharmProjects/KinaseModulation/CosmicFiles/cosmic_filtered_all_psites.csv'
cosmic_percent_rankCSV = 'C:/Users/eelias13/PycharmProjects/KinaseModulation/CosmicFiles/cosmic_percent_rank.csv'
#cosmic_four_samplesCSV = 'C:/Users/eelias13/PycharmProjects/KinaseModulation/CosmicFiles/cosmic_sample.csv'

gtex_metaCSV = '/home/eelias13/KinaseModulation/gtexFiles/GTEX_metadata.csv'
gtex_tpmTSV = '/home/eelias13/KinaseModulation/gtexFiles/GTEx_v1.1.9_Median_tpm_by_SMTS.tsv'
tissues_mappedCSV = "/home/eelias13/KinaseModulation/gtexFiles/cosmic_tissue_types_mapped.csv"
cosmic_tissuesCSV = '/home/eelias13/KinaseModulation/gtexFiles/cosmic_tissues.csv'
tpm_filteredCSV = 'C/home/eelias13/KinaseModulation/gtexFiles/tpm_filtered.csv'
sty_kinase_tpmCSV = '/home/eelias13/KinaseModulation/gtexFiles/sty_kinase_tpm_mapped.csv'
cosmic_tpmCSV = '/home/eelias13/KinaseModulation/gtexFiles/cosmic_variant_tpm_mapped.csv'
poiCSV = '/home/eelias13/KinaseModulation/gtexFiles/proteins_of_interest.csv'
psite_tpmCSV = '/home/eelias13/KinaseModulation/gtexFiles/psite_tpm_mapped.csv'

ser_thr_kinasesTSV = '/home/eelias13/KinaseModulation/kinaseFiles/mapped_ser_thr_kinases.tsv'
tyr_kinasesTSV = '/home/eelias13/KinaseModulation/kinaseFiles/mapped_tyr_kinases.tsv'
ser_thr_tyr_kinasesCSV = '/home/eelias13/KinaseModulation/kinaseFiles/mapped_ser_thr_tyr_kinases.csv'

HT_psitesTSV = '/home/eelias13/KinaseModulation/final_phosphosites_mapped_to_canonical_proteome-merged.tsv'

# data frame variables
# proteome, kinase, and phosphosite dfs
proteome = readCSV(proteomeCSV)
final_phosphosites = readCSV(final_phosphositesCSV)
#substrate_yaffe = readCSV(substrate_yaffeCSV)
#substrate_psp = readCSV(substrate_pspCSV)
all_psites_mapped_ensembl = readCSV(all_psites_mapped_ensemblCSV)

# pssm df
st_pssm = readCSV(pssmCSV)
y_pssm = readCSV(tyrPSSMCSV)

# cosmic dfs (unaltered)
cosmic_full = readTSV(cosmicTSV)
cosmic_head = readCSV(cosmic_headCSV) #tester df (first 1000 lines)

# cosmic altered dfs
cosmic_all_psites = readCSV(cosmic_all_psitesCSV)
cosmic_phos_muts = readCSV(cosmic_phos_mutsCSV)
cosmic_scores = readCSV(cosmic_scoresCSV)
cosmic_percent_rank = readCSV(cosmic_percent_rankCSV)
# cosmic_samples_dict = {
#     "#Gene": ["SGO1", "GCSAML", "SRC", "NACA"],
#     "UniProt_accession": ["Q5FBB7", "Q5JQS6", "P12931", "E9PAV3"],
#     "UniProt_position": [448, 116, 439, 1413],
#     "residue_WT": ["S", "S", "Y", "K"], 
#     "residue_variant": ["F", "Y", "S", "Q"]
# }
# cosmic_samples = pd.DataFrame(cosmic_samples_dict)
# PM_cosmic_sample = readCSV('PM_cosmic_sample.csv')
# sample_cosmic_scores = readCSV('cosmic_sample_kinase_scores.csv')

# phosphosites df
proteome = readCSV(proteomeCSV)

# tissues
gtex_meta = readCSV(gtex_metaCSV)
gtex_tpm = readTSV(gtex_tpmTSV)
tissues_mapped = readCSV(tissues_mappedCSV)
cosmic_tissues = readCSV(cosmic_tissuesCSV)
tpm_filtered = readCSV(tpm_filteredCSV)
sty_kinase_tpm = readCSV(sty_kinase_tpmCSV)
cosmic_tpm = readCSV(cosmic_tpmCSV)
poi = readCSV(poiCSV)
psite_tpm = readCSV(psite_tpmCSV)

# kinases mapped
ser_thr_kinases = readTSV(ser_thr_kinasesTSV)
tyr_kinases = readTSV(tyr_kinasesTSV)
ser_thr_tyr_kinases = readCSV(ser_thr_tyr_kinasesCSV)

# HT psites
ht_psites = readTSV(HT_psitesTSV)

# adjust data types and indices in df as needed
# command example: df.col = df.col.astype('str')


# remove all non-digits from phosphosite column
def removeDigits(df):
    df["PHOSPHOSITE"] = df["PHOSPHOSITE"].apply(lambda x: re.sub(r'\D', "", x))
    return df


# method to return the uniprot entry ID given a substrate sequence
# provide the given sequence and specify data frame it will be compared to
def getUniprotID(sequence, df):
    row = df.index.get_loc(df[df['SITE_+/-7_AA'] == sequence].index[0])
    access_id = df.at[row, 'ACC_ID']
    return access_id


# method to return the phosphorylation site of a given sequence
def getPhosphosite(sequence, df):
    row = df.index.get_loc(df[df["SITE_+/-7_AA"] == sequence].index[0])
    site = df.loc[row, "PHOSPHOSITE"]
    return site


# method to return the full protein sequence from the proteome dataframe given an accession id
def getSeq(uniprot_id):
    uniprot_id = str(uniprot_id)
    if (proteome == uniprot_id).any().any():
        row = proteome.index.get_loc(proteome[proteome['ACC_ID'] == uniprot_id].index[0])
        sequence = proteome.at[row, 'PROTEIN_SEQUENCE']
        return sequence
    else:
        return None


# method to check a given sequence against the proteome
# returns None if the site is NOT valid
# if phosphorylation site in the full proteome contains ser/thr, return the sequence
# accepts a target sequence and df target sequence was obtained from
def checkProteome(sequence, df, site_df):
    currentID = getUniprotID(sequence, df)
    site = getPhosphosite(sequence, site_df)
    # get residue at the equivalent location in the full proteome
    row = proteome.index.get_loc(proteome[proteome["ACC_ID"] == currentID].index[0])
    proteome_sequence = proteome.at[row, "PROTEIN_SEQUENCE"]
    # compare residues
    compareTo = proteome_sequence[site - 1]
    if compareTo == "S" or compareTo == "T":
        return sequence
    else:
        return None


# method to generate a dictionary with a substrate and a list of all sequences
# accepts a dataframe of phosphosites
def substrateSequences(phosphosite_df):
    seq_ID = {}
    for i in range(len(phosphosite_df)):
        uniprot_id = phosphosite_df.loc[i, "ACC_ID"]
        sequence = phosphosite_df.loc[i, "SITE_+/-7_AA"]
        if uniprot_id not in seq_ID:
            seq_ID[uniprot_id] = [sequence]
        else:
            seq_ID[uniprot_id].append(sequence)
    return seq_ID


# function to convert a dataframe into a csv
def to_csv(df, filename):
    df.to_csv(f'{filename}.csv')
    return


# function to check if there is a phosphorylatable residue within +/-5 AA of a residue
# checks if there is a S/T within 4 AA before the residue or 5 AA after the residue
# returns boolean value
def find_st(position, uniprot_id):
    sequence = str((getSeq(uniprot_id))).upper()
    upstream = sequence[(position - 4):position]
    downstream = sequence[position:(position + 5)]
    if sequence[position - 1] == 'S' or sequence[position - 1] == 'T':
        return 'phosphosite'
    elif 'S' in upstream or 'T' in upstream:
        return 'upstream'
    elif 'S' in downstream or 'T' in downstream:
        return 'downstream'
    else:
        return None


# function to check if there is a phosphorylatable residue within +/-5 AA of a residue
# checks if there is a S/T within 4 AA before the residue or 5 AA after the residue
# returns boolean value
def find_y(position, uniprot_id):
    sequence = str((getSeq(uniprot_id))).upper()
    upstream = sequence[(position - 4):position]
    downstream = sequence[position:(position + 5)]
    if sequence[position - 1] == 'Y':
        return 'phosphosite'
    elif 'Y' in upstream:
        return 'upstream'
    elif 'Y' in downstream:
        return 'downstream'
    else:
        return None


