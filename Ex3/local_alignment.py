import numpy as np

from tests.matrices import MATRICES

class LocalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        self.string1 = string1
        self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        print(string1)
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int)
        self.align()

    def align(self):
        """
        Align given strings using the Smith-Waterman algorithm.
        NB: score matrix and the substitution matrix are different matrices!
        """
        print("sdfsfdfsdfdsfsf")
        scr_mat = self.score_matrix
        
        gap = self.gap_penalty
        #gap = -2
        str1 = self.string1
        str2 = self.string2
        scr_mat[0][:] = 0
        scr_mat[:][0] = 0

        ref_dict = {}
        sub_mat = self.substitution_matrix
        for j in range(1, len(str2) + 1 ):
            for i in range(1, len(str1) + 1):
                sub_val = sub_mat[str2[j-1]][str1[i-1]]
                hor = scr_mat[j][i-1] + gap 
                ver = scr_mat[j-1][i] + gap
                if str1[i-1] == str2[j-1]:
                    diag = scr_mat[j-1][i-1] +  sub_val
                else:
                    #diag = scr_mat[j-1][i-1] + mismatch+ sub_val
                    diag = scr_mat[j-1][i-1] + sub_val
                val_list = [hor, ver, diag]
                m = max(val_list)
                max_pos = [i for i, j in enumerate(val_list) if j == m]
                tup = (j, i)
                ref_dict[tup] = max_pos
                m = max(0, m)
                scr_mat[j][i] = m
        #max_sc_ver = np.argwhere(scr_mat[: ,:] == np.amax(scr_mat[:, :])).flatten().tolist()
        max_sc_ver = np.argwhere(scr_mat[: ,:] == np.amax(scr_mat[:, :])).tolist()
        #max_sc_hor = np.argwhere(scr_mat[-1, :] == np.amax(scr_mat[-1, :])).flatten().tolist()
        """
        max_sc_ver = np.argwhere(scr_mat == np.amax(scr_mat)).flatten().tolist()   
        if not any(isinstance(i, list) for i in max_sc_ver):
            max_sc_ver = [max_sc_ver]
        """
        len_str1 = len(str1)
        len_str2 = len(str2)
        seqs = []
        str_indx = {'str1':[], 'str2': []}
        for indx in max_sc_ver:
        #for indx in max_sc_hor:
            print(indx)
            isDone = False
            while not isDone and indx[0] != 0 and indx[1] != 0:
                count = 0
                seq = []
                main_seq = []
                i = 0
                for i in range(len_str1):
                    if len(seq) == 0:
                        #pos = ref_dict[(len_str2 -1, indx)][0]
                        pos = ref_dict[(indx[0], indx[1])][0]
                        #pos = ref_dict[(indx[0], indx[1])][0]
                        if pos == 2:
                            #seq.append(str2[indx[0] - 1])
                           
                            seq.append(str2[indx[0]-1])
                            str_indx['str2'].append(indx[0]-1)
                            #seq.append(str2[indx- 1])
                            #p1 = indx[0] -1 
                            p1 = indx[0] - 1
                            #p2 = len_str1 -1 
                            p2 = indx[1] - 1
                        elif pos == 0:
                            seq.append('-')
                            p1 = indx[0]
                            p2 = indx[1] - 1
                            #p1 = len_str2  
                            #p2 = len_str1 - 1 

                        elif pos == 1:
                            #p1 = len_str2 - 1 
                            #p2 = len_str1
                            p1 = indx[0] - 1
                            p2 = indx[1] 
                            seq.append('-')
                        print("sdfsfdfsdfdsfsf")
                        main_seq.append(str1[indx[1] - 1])
                        str_indx['str1'].append(indx[1] - 1)
                    elif scr_mat[p1][p2] != 0:
                    #if(True):
                        pos = ref_dict[(p1, p2)]
                        if len(pos) > 1:
                            count += 1
                            pos = pos[0]
                            ref_dict[(p1, p2)].remove(pos)
                        else:
                            pos = pos[0]
                        if pos == 2:
                            seq.append(str2[p1  - 1])
                            str_indx['str2'].append(p1  - 1)
                            p1 = p1 - 1
                            p2 = p2 - 1

                        elif pos == 0:
                            p1 = p1 
                            p2 = p2 - 1
                            seq.append('-')
                        elif pos == 1:
                            p1 = p1 
                            p2 = p2 - 1
                            seq.append('-')
                        main_seq.append(str1[p2])
                        str_indx['str1'].append(p2)
                        
                if count > 0:
                    isDone = False
                else:
                    isDone = True          

                seq.reverse()
                seqs.append(seq)
                main_seq.reverse()
        
       # print("sdfsfdfsdfdsfsf")
       # if len(seqs) > 1:
        tot_scores = {}
        sub_mat = self.substitution_matrix
        for seq in seqs:
            tot_score = 0
            for k in range(len(seq)):
                if seq[k] == '-':
                    tot_score += self.gap_penalty
                else:
                    tot_score += sub_mat[str1[k]][seq[k]]   
            tot_scores[''.join(seq)]  = tot_score   
        if len(seq) > 2:
            max_value = max(tot_scores.values())
            self.best_score = max_value
            all_seqs = [k for k,v in tot_scores.items() if v == max_value]
            final_seqs = []
            for final in all_seqs:
                final_seqs.append((''.join(main_seq), final))
            self.alignments = final_seqs[0]
            self.str_indx = str_indx
        else:
            self.alignments = ("", "")
            self.str_indx = {}
       # else:
        #    final_seqs  = [(str1, ''.join(seqs[0]))]
        return self.alignments

    def has_alignment(self):
        """
        :return: True if a local alignment has been found, False otherwise
        """
        if self.alignments[0] != '':
            return True
        else:
            return False

    def get_alignment(self):
        """
        :return: alignment represented as a tuple of aligned strings
        """
        return self.alignments
        #return ('DAC', 'DAC')

    def is_residue_aligned(self, string_number, residue_index):
        """
        :param string_number: number of the string (1 for string1, 2 for string2) to check
        :param residue_index: index of the residue to check
        :return: True if the residue with a given index in a given string has been alined
                 False otherwise
        """
        res = False
        if string_number == 1:
            if residue_index in self.str_indx['str1']:
                res = True
        elif string_number == 2:
            if residue_index in self.str_indx['str2']:
                res = True     
        return res


"""
input = {"strings": ["ARNDCEQGHI", "DDCEQHG"],
        "gap_penalty": -6,
        "matrix": "blosum",
        "alignment": ["NDCEQGH", "DDCEQ-H"],
        "residue_aligned_on_first": [[1, 0, False], [1, 3, True]],
        "residue_aligned_on_second": [[2, 0, True], [2, 6, False]]}
g =  LocalAlignment(input["strings"][0], input["strings"][1], input["gap_penalty"], MATRICES[input['matrix']])
g.is_residue_aligned(1, 0)
"""