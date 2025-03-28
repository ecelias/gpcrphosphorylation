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
proteomeCSV = "/Users/elizabethelias/StJude/gpcrphosphorylation/datafiles/proteome.csv"
pssmCSV = "/Users/elizabethelias/StJude/gpcrphosphorylation/datafiles/ser_thr_pssm.csv"
tyrPSSMCSV = '/Users/elizabethelias/StJude/gpcrphosphorylation/datafiles/tyr_pssm.csv'
final_phosphositesCSV = "/Users/elizabethelias/StJude/gpcrphosphorylation/datafiles/final_mapped_phosphosites.csv"


# data frame variables
# proteome, kinase, and phosphosite dfs
proteome = readCSV(proteomeCSV)
final_phosphosites = readCSV(final_phosphositesCSV)

# pssm df
st_pssm = readCSV(pssmCSV)
y_pssm = readCSV(tyrPSSMCSV)

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


