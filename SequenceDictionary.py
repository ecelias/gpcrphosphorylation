import Aim1

substrate_yaffe = Aim1.substrate_yaffe
substrate_psp = Aim1.substrate_psp

# generate a dictionary where the key is ACC_ID and value is the full substrate sequence
# for yaffe dataset
def yaffe_acc_dict():
    yaffe_acc_dict = {}
    for i, row in substrate_yaffe.iterrows():
        row = i
        uniprot_id = substrate_yaffe.loc[row, "ACC_ID"]
        substrate_seq = Aim1.getSeq(uniprot_id)
        if substrate_seq is not None:
            yaffe_acc_dict[uniprot_id] = substrate_seq
        else:
            print(f"Warning: Sequence not found for UniProt ID {uniprot_id}")
    return yaffe_acc_dict


# for psp datset
def psp_acc_dict():
    psp_acc_dict = {}
    for i, row in substrate_psp.iterrows():
        row = i
        uniprot_id = substrate_psp.loc[row, "ACC_ID"]
        substrate_seq = Aim1.getSeq(uniprot_id)
        if substrate_seq is not None:
            psp_acc_dict[uniprot_id] = substrate_seq
        else:
            print(f"Warning: Sequence not found for UniProt ID {uniprot_id}")
    return psp_acc_dict