import Aim1
import pickle

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

#mapped_ht = all_psites_dict(Aim1.ht_psites)
#print(mapped_ht)
#pickle.dump(mapped_ht, open("pickle_ht_psites_dict.p", 'wb'))
mapped_ht_psites = pickle.load(open("pickle_ht_psites_dict.p", "rb"))


