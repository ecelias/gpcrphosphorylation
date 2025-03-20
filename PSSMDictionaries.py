import Aim1
import pickle

st_pssm = Aim1.st_pssm
y_pssm = Aim1.y_pssm


# create a dictionary where k = (gene, position, residue) and v = PSSM score
def pssm_dict(df):
    dict = {}
    for index, row in df.iterrows():
        gene_name = row["GENE"]
        for column, item in df.items():
            residue_site = str(column)
            residue = residue_site[-1]
            site = residue_site[:-1]
            score = item[index]
            if not residue_site == 'GENE':
                k = gene_name + ',' + site + ',' + residue
                k = str(k)
                dict[k] = score
    return dict


#st_pssm_dictionary = pssm_dict(st_pssm)
#pickle.dump(st_pssm_dictionary, open("pickle_st_pssm_dictionary.p", 'wb'))
st_pssm_dict = pickle.load(open("pickle_st_pssm_dictionary.p", "rb"))

#y_pssm_dictionary = pssm_dict(y_pssm)
#pickle.dump(y_pssm_dictionary, open("pickle_y_pssm_dictionary.p", 'wb'))
y_pssm_dict = pickle.load(open("pickle_y_pssm_dictionary.p", "rb"))

