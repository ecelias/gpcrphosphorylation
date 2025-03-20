import pandas as pd 
import pickle
import math
from scipy.stats import norm

class getRank:
    psite_marker = "*"

    # pssm data
    pssmCSV = "pssmFiles/pssm.csv"
    tyrPSSMCSV = 'pssmFiles/tyrosine_pssm.csv'

    def __init__(self, sequence, score_csv=False, pr_csv=False, use_ht=False):

        """ 
        Parameters: 
            sequence = amino acid sequence of any length which phosphosite residues followed by a *
            score_csv = (optional) save a CSV file containing the PSSM scores for each kinase. Default is False.
            pr_csv = (optional) save a CSV file containing the percent rank for each kinase. Default is False.
            use_ht = (optional) uses the high throughput phosphosites to calculate the percent rank. Default is False

        Example: test = getRank(sequence="LEPLKS*LDQELEPLKY*LDQE", score_csv=True, pr_csv=True, use_ht=True) """

        self.sequence = sequence
        self.sequence_list = []

        # Load the pssm data files
        self.st_pssm = pd.read_csv(self.pssmCSV)
        self.y_pssm = pd.read_csv(self.tyrPSSMCSV)
        
        # Extract kinase lists
        self.st_kinase_list = self.st_pssm['GENE'].tolist()
        self.y_kinase_list = self.y_pssm['GENE'].tolist()
        self.all_kinases = self.st_kinase_list + self.y_kinase_list

        # Load pickled pssm dictionaries
        try:
            self.st_pssm_dict = pickle.load(open("pickle_st_pssm_dictionary.p", "rb"))
            self.y_pssm_dict = pickle.load(open("pickle_y_pssm_dictionary.p", "rb"))
        except FileNotFoundError:
            print("Pickle files not found.")
        except Exception as e:
            print(f"Error loading pickle files: {e}")

        # Load score distributions for percent rank
        self.use_ht = use_ht
        self.sd_st = pd.read_csv('st_stats_df.csv', index_col=0)
        self.sd_y = pd.read_csv('y_stats_df.csv', index_col=0)
        self.ht_sd_st = pd.read_csv('ht_st_stats_df.csv', index_col=0)
        self.ht_sd_y = pd.read_csv('ht_y_stats_df.csv', index_col=0)

        # boolean to indicate if you want score and percent rank df, default is false
        self.score_csv = score_csv
        self.pr_csv = pr_csv
    
    def getSeq(self):
        print(self.sequence)
        return self.sequence

    def getSaveScores(self):
        return self.score_csv

    def getSavePR(self):
        return self.pr_csv

    def useHT(self):
        return self.use_ht
    
    def getKinaseLists(self):
        return self.all_kinases, self.st_kinase_list, self.y_kinase_list
    
    def getPSSMDict(self):
        return self.st_pssm_dict, self.y_pssm_dict

    def getScoreDistributions(self):
        ht = self.useHT()
        if ht is False:
            return self.sd_st, self.sd_y
        else:
            return self.ht_sd_st, self.ht_sd_y

    def getSeqList(self):
        return self.sequence_list

    def updateSeqList(self, sequence):
        self.sequence_list.append(sequence)

    def countPSites(self):

        """ Counts the number of phosphosites that exist in the sequence used to create the 
        getRank object. Returns both the number of sequences and a list of the indices. The 
        list of indices contains the index of all phosphosites in a sequence if the 
        phosphosite marker was removed from the string. """

        seq = self.getSeq()
        psiteIdxs = []
        numSites = 0
        if seq.find(self.psite_marker) < 0:
            raise Exception(f"No phosphosites marked. Please mark psites with {self.psite_marker}")
        else:
            numSites = seq.count(self.psite_marker)
            j = 1
            for i in range(len(seq)):
                if seq[i] == self.psite_marker:
                    psiteIdxs.append(i - j)
                    j += 1
        return numSites, psiteIdxs

    def extractFlanking(self):

        """ Uses the countPsites function to extract the flanking region residues surrounding
        each phosphosite marked in the sequence. If necessary, adjusts the length of the subsequence containing
        the flanking region to have a length of 10 by adding underscores to the front or back of the
        sequence depending on the location of the phosphosite. """

        self.countPSites()

        seq = self.getSeq()
        seq = seq.replace(self.psite_marker, '')
        seq_length = len(seq)
        num_sites, site_idxs  = self.countPSites()
        for site in site_idxs:
            if site < 5:
                currentSeq = seq[:(site+5)]
                if len(currentSeq) < 10:
                    currentSeq.rjust(10, '_')
            elif 5 <= site <= (seq_length - 4):
                currentSeq = seq[(site - 5):(site+5)]
            elif site > (seq_length - 4):
                currentSeq = seq[(site - 5):]
                if len(currentSeq) < 10:
                    currentSeq.ljust(10, '_')
            self.updateSeqList(currentSeq)
    
    def which_family(self, seq):

        """Helper function that determines the family of kinases which phosphorylate
        a sequence. Returns 1 if Ser/Thr family kinase, returns 2 if Tyr family kinase, and returns
        -1 if the phosphosite doesn't exist."""

        psite = seq[5].lower()
        return {"s": 1, "t": 1, "y": 2}.get(psite, -1)
        
    # helper function that accepts a kinase, site, and residue and returns the PSSM
    def get_pssm_st(self, kinase, site, residue):

        """ Helper function that returns the PSSM score for a residue for Serine/Threonine family kinases. Must
        have the correct path to the pickle_st_pssm_dictionary.p file."""

        st_pssm_dict = self.getPSSMDict()[0]

        key = str(kinase) + ',' + str(site) + ',' + str(residue)
        score = st_pssm_dict.get(key)
        if score is not None:
            return float(score)


    def get_pssm_y(self, kinase, site, residue):
        
        """ Helper function that returns the PSSM score for a residue for Tyrosine family kinases. Must
        have the correct path to the pickle_y_pssm_dictionary.p file."""

        y_pssm_dict = self.getPSSMDict()[1]

        key = str(kinase) + ',' + str(site) + ',' + str(residue)
        score = y_pssm_dict.get(key)
        if score is not None:
            return float(score)
    
    def score_st(self, sequence):

        """Returns the PSSM score for an amino acid sequence if the phosphosite is a serine or threonine for all
        Ser/Thr family kinases.

        The total PSSM is obtained by finding the product of the PSSM value for each residue in the sequence and
        taking the log2 value of that score. Total PSSM score will always equal 0 for Tyr family kinases. """

        all_kinases = self.getKinaseLists()[0]
        st_kinase_list = self.getKinaseLists()[1]
        y_kinase_list = self.getKinaseLists()[2]

        st_score_dict = {}
        sequence = str(sequence)

        for kinase in all_kinases:
            st_current_site = -5
            st_seq_index = 0
            st_product = 1
            while st_current_site < 5:
                if st_seq_index > len(sequence) - 1:
                    break
                if kinase in y_kinase_list:
                    st_product = 1
                elif kinase in st_kinase_list and st_current_site != 0:
                    st_current_residue = sequence[st_seq_index]
                    if st_current_residue == "_":
                        st_product *= 1
                    elif self.get_pssm_st(kinase, st_current_site, st_current_residue) is None:
                        st_product *= 1
                        print(f"Warning: Residue {st_current_residue} at {st_current_site} does not have PSSM score")
                    else:
                        st_product *= self.get_pssm_st(kinase, st_current_site, st_current_residue)
                st_seq_index += 1
                st_current_site += 1
            st_score_dict[kinase] = math.log2(st_product)
        return st_score_dict

    def score_y(self, sequence):

        """Returns the PSSM score for an amino acid sequence if the phosphosite is a tyrosine for all
        Tyr family kinases.

        The total PSSM is obtained by finding the product of the PSSM value for each residue in the sequence and
        taking the log2 value of that score. Total PSSM score will always equal 0 for Ser/Thr family kinases. """

        all_kinases = self.getKinaseLists()[0]
        st_kinase_list = self.getKinaseLists()[1]
        y_kinase_list = self.getKinaseLists()[2]        

        y_score_dict = {}
        sequence = str(sequence)

        for kinase in all_kinases:
            y_current_site = -5
            y_seq_index = 0
            y_product = 1
            while y_current_site < 5:
                if y_seq_index > len(sequence) - 1:
                    break
                if kinase in st_kinase_list:
                    y_product = 1
                elif kinase in y_kinase_list and y_current_site != 0:
                    y_current_residue = sequence[y_seq_index]
                    if y_current_residue == "_":
                        y_product *= 1
                    elif self.get_pssm_y(kinase, y_current_site, y_current_residue) is None:
                        y_product *= 1
                        print(f"Warning: Residue {y_current_residue} at {y_current_site} does not have PSSM score")
                    else:
                        y_product *= self.get_pssm_y(kinase, y_current_site, y_current_residue)
                y_seq_index += 1
                y_current_site += 1
            y_score_dict[kinase] = math.log2(y_product)
        return y_score_dict        

    def score_seqs(self):

        """Returns a dataframe that contains the total PSSM score for all kinases."""

        self.extractFlanking()

        all_kinases = self.getKinaseLists()[0]
        st_kinase_list = self.getKinaseLists()[1]
        y_kinase_list = self.getKinaseLists()[2]

        column_names = ["SEQ", "PSITE_RESIDUE"] + all_kinases
        score_df = pd.DataFrame(columns=column_names)
        seq_list = self.getSeqList()
        print(seq_list)

        for sequence in seq_list:
            family = self.which_family(sequence)
            psite_res = None
            score_dict = {}
            if family is None or family == -1:
                raise Exception("No phosphosite found. Cannot score sequence.")
            if family == 1:
                # score st
                score_dict = self.score_st(sequence)
                psite_res = "ST"
            elif family == 2:
                # score y 
                score_dict = self.score_y(sequence)
                psite_res = "Y"
            info = {"SEQ": sequence, "PSITE_RESIDUE": psite_res}
            new_column = {**info, **score_dict}
            score_df = score_df.append(new_column, ignore_index=True)

        if self.getSaveScores() is True:
            score_df.to_csv("kinase_scores.csv", index=False)     

        return score_df

    def percent_rank(self):

        """ Returns a dataframe containing the percent rank the likelihood of each kinase to phosphorylate each sequence.
        If the kinase belongs to a family that would be unable to phosphorylate a sequence, then the percent rank is 
        set to 0. """

        all_kinases = self.getKinaseLists()[0]
        st_kinase_list = self.getKinaseLists()[1]
        y_kinase_list = self.getKinaseLists()[2]

        sd_st = self.getScoreDistributions()[0]
        sd_y = self.getScoreDistributions()[1]
        score_df = self.score_seqs()
        percent_dict = {}

        for index, row in score_df.iterrows():
            kinase_dict = {}
            psite_res = str(row['PSITE_RESIDUE'])
            kinase_dict["PSITE_RESIDUE"] = psite_res

            for kinase in all_kinases:
                if kinase in st_kinase_list and psite_res == "ST":
                    mean = sd_st.loc["mean", kinase]
                    std = sd_st.loc["std", kinase]
                    score = row[kinase]
                    z_score = (score - mean) / std
                    percent = norm.cdf(z_score) * 100
                    kinase_dict[kinase] = percent
                elif kinase in y_kinase_list and psite_res == "Y":
                    mean = sd_y.loc["mean", kinase]
                    std = sd_y.loc["std", kinase]
                    score = row[kinase]
                    z_score = (score - mean) / std
                    percent = norm.cdf(z_score) * 100
                    kinase_dict[kinase] = percent
                else:
                    kinase_dict[kinase] = 0

            percent_dict[str(row['SEQ'])] = kinase_dict

        percent_df = pd.DataFrame.from_dict(percent_dict, orient='index')
        percent_df = percent_df.reset_index() 
        pd.options.display.float_format = '{:.3f}'.format   

        if self.getSavePR() is True:
            percent_df.to_csv("percent_rank.csv", index=False) 

        return percent_df


#test = getRank("LEPLKS*LDQELEPLKY*LDQE", True, True)
