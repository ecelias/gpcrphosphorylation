import Aim1
import pickle
import pandas as pd

final_sites = Aim1.final_phosphosites

# create a dictionary where the key is [uniprot id, phosphorylated residue] and the value is [+/- 5 AA]
# all mapped phosphosites dataset
def all_psites_dict(psites_df):
    all_psites_dict = {}
    for index, row in psites_df.iterrows():
        current_gene = row["#UniProt_accession"]
        print(current_gene)
        current_gene = str(current_gene)
        if Aim1.getSeq(current_gene) is not None:
            site = row["position"]
            residue = row['residue']
            k = str(current_gene) + ',' + str(site)
            current = Aim1.getSeq(current_gene)
            residues = ''
            if site - 6 < 0:
                empty = abs(site - 6)
                remaining = 6 - empty
                i = 0
                while i < empty:
                    residues = residues + '_'
                    i += 1
                residues = residues + current[(site - remaining):(site + 4)]
            elif site + 4 > len(current):
                remaining = (site + 4) - len(current)
                if remaining > 10:
                    residues = '__________'
                else:
                    empty = 4 - remaining
                    residues = current[(site - 6):(site + remaining)]
                    i = 0
                    while i < empty:
                        if len(residues) < 11:
                            residues = residues + '_'
                            i += 1
                        else:
                            break
            else:
                residues = current[(site - 6):(site + 4)]
            all_psites_dict[k] = residues + ',' + str(residue)
        else:
            print(f"Warning: No phosphosite found for UniProt ID {current_gene}")
    return all_psites_dict


#all_mapped_dict = all_psites_dict()
#pickle.dump(all_mapped_dict, open("pickle_all_psites_dict.p", 'wb'))
all_mapped_psites = pickle.load(open("pickle_all_psites_dict.p", "rb"))

def to_df():
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(all_mapped_psites, orient='index', columns=['Sequence_Info'])

    # Split the index into Protein ID and Position
    df.index.name = 'Protein_Position'
    df = df.reset_index()
    df[['UniProtID', 'Position']] = df['Protein_Position'].str.split(',', expand=True)
    df[['Sequence', 'Phosphosite']] = df['Sequence_Info'].str.split(',', expand=True)

    # Reorder columns and drop the intermediate
    df = df[['UniProtID', 'Position', 'Sequence', 'Phosphosite']]
    # Ensure Position is an integer
    df['Position'] = df['Position'].astype(int)

    # Sort by both Protein and Position
    mapped_psites = df.sort_values(by=['UniProtID', 'Position']).reset_index(drop=True)

    # Display the DataFrame
    print(mapped_psites)
    return mapped_psites

mapped_psites_df = to_df()
mapped_psites_df['Position'] = mapped_psites_df['Position'].astype(int)

def psites_in_flanking():
    mapped_psites_df['FlankPSiteCount'] = 0
    mapped_psites_df['FlankPSiteIndices'] = [[] for _ in range(len(mapped_psites_df))] 

    # Loop through each row
    for idx, row in mapped_psites_df.iterrows():
        uid = row['UniProtID']
        center_pos = row['Position']
        seq = row['Sequence']
        phosphosite_res = row['Phosphosite']

        # The central phosphosite is at index 5 in the sequence
        # So, the protein position corresponding to seq[0] is:
        start_protein_pos = center_pos - 5

        # Get all other canonical phosphosites in same protein (excluding self)
        same_protein = mapped_psites_df[(mapped_psites_df['UniProtID'] == uid) & (mapped_psites_df.index != idx)]

        # Prepare to collect valid phosphosite indices in the flanking region
        flank_indices = []
        flank_indices_str = ""

        for _, other_row in same_protein.iterrows():
            other_pos = other_row['Position']
            offset = other_pos - start_protein_pos  # map protein position to sequence index

            # Check if within the flanking region (0–4 before, 6–9 after)
            if 0 <= offset < len(seq) and offset != 5 and (offset <= 4 or offset >= 6):
                flank_indices.append(offset)
                if flank_indices_str == "":
                    flank_indices_str = str(offset)
                else:
                    flank_indices_str = flank_indices_str + ";" + str(offset)
    
        mapped_psites_df.at[idx, 'FlankPSiteCount'] = len(flank_indices)
        mapped_psites_df.at[idx, 'FlankPSiteIndices'] = flank_indices_str


    return mapped_psites_df

#mapped_psites_with_flanking_info = psites_in_flanking()
#Aim1.to_csv(mapped_psites_with_flanking_info, "mapped_psites_with_flanking_info")