import pandas as pd

# define path to gunzipped data
file_path = 'AlphaMissense_aa_substitutions.tsv'
# define a list of accession numbers of proteins of interest
acc = ['P21728', 'P14416', 'P35462', 'P21917', 'P21918']
# read in 6 GB AM dataset (This takes a while)
AM = pd.read_csv(file_path, sep='\t', skiprows=3)
# reset the index to the accession number
AM.set_index('uniprot_id', inplace=True)
# write a separate .tsv file for all accession numbers of a proteins of interest
for i in acc:
    AM.loc[i].to_csv(i + '.tsv', sep='\t', index=False)