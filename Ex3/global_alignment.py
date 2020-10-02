import numpy as np
from tests.matrices import MATRICES


class GlobalAlignment:
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
        self.substituion_matrix = matrix
        #print(matrix)
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int)
        self.align()

    def align(self):
        """
        Align given strings using the Needleman-Wunsch algorithm,
        store the alignments and the score matrix used to compute those alignments.
        NB: score matrix and the substitution matrix are different matrices!
        """

        scr_mat = self.score_matrix
        match = 1
        mismatch = -1
        gap = self.gap_penalty
        #gap = -2
        str1 = self.string1
        str2 = self.string2
        for i in range(len(str1) + 1):
            if i == 0:
                scr_mat[0][i] = 0
            else:

                scr_mat[0][i] = scr_mat[0][i-1] + gap

        for i in range(len(str2) + 1):
            if i == 0:
                scr_mat[i][0] = 0
            else:

                scr_mat[i][0] = scr_mat[i-1][0] + gap
        ref_dict = {}
        sub_mat = self.substituion_matrix
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
                scr_mat[j][i] = m
        max_sc_ver = np.argwhere(scr_mat[: ,-1] == np.amax(scr_mat[:, -1])).flatten().tolist()
        max_sc_hor = np.argwhere(scr_mat[-1, :] == np.amax(scr_mat[-1, :])).flatten().tolist()
        """
        max_sc_ver = np.argwhere(scr_mat == np.amax(scr_mat)).flatten().tolist()   
        if not any(isinstance(i, list) for i in max_sc_ver):
            max_sc_ver = [max_sc_ver]
        """
        len_str1 = len(str1)
        len_str2 = len(str2)
        seqs = []
        for indx in max_sc_ver:
        #for indx in max_sc_hor:
            
            isDone = False
            while not isDone:
                count = 0
                seq = []
                i = 0
                for i in range(len_str1):
                    
                    if len(seq) == 0:
                        #pos = ref_dict[(len_str2 -1, indx)][0]
                        pos = ref_dict[(indx, len_str1)][0]
                        #pos = ref_dict[(indx[0], indx[1])][0]
                        if pos == 2:
                            #seq.append(str2[indx[0] - 1])
                            """
                            if len_str1 - 1 > len_str2 - 1:
                                seq.append("-")
                            else:
                                seq.append(str2[len_str1- 1])
                            """
                            seq.append(str2[indx-1])
                            #seq.append(str2[indx- 1])
                            p1 = len_str2 -1 
                            #p1 = indx[0] - 1
                            p2 = len_str1 -1 
                            #p2 = indx[1] - 1
                        elif pos == 0:
                            seq.append('-')
                            #p1 = indx[0]
                            #p2 = indx[1] - 1
                            p1 = len_str2  
                            p2 = len_str1 - 1 

                        elif pos == 1:
                            p1 = len_str2 - 1 
                            p2 = len_str1
                            #p1 = indx[0] - 1
                            #p2 = indx[1] 
                            seq.append('-')
                    else:
                        pos = ref_dict[(p1, p2)]
                        if len(pos) > 1:
                            count += 1
                            pos = pos[0]
                            ref_dict[(p1, p2)].remove(pos)
                        else:
                            pos = pos[0]
                        if pos == 2:
                            seq.append(str2[p1  - 1])
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
                        
                if count > 0:
                    isDone = False
                else:
                    isDone = True          

                seq.reverse()
                seqs.append(seq)


       # if len(seqs) > 1:
        tot_scores = {}
        sub_mat = self.substituion_matrix
        for seq in seqs:
            tot_score = 0
            for i in range(len_str1):
                if seq[i] == '-':
                    tot_score += self.gap_penalty
                else:
                    tot_score += sub_mat[str1[i]][seq[i]]   
            tot_scores[''.join(seq)]  = tot_score   

        max_value = max(tot_scores.values())
        self.best_score = max_value
        all_seqs = [k for k,v in tot_scores.items() if v == max_value]
        final_seqs = []
        for final in all_seqs:
            final_seqs.append((str1, final))
        self.alignments = final_seqs
       # else:
        #    final_seqs  = [(str1, ''.join(seqs[0]))]
        
        return final_seqs



    def get_best_score(self):
        """
        :return: the highest score for the aligned strings, int

        """
        return self.best_score

    def get_number_of_alignments(self):
        """
        :return: number of found alignments with the best score
        """
        return len(self.alignments)

    def get_alignments(self):
        """
        :return: list of alignments, where each alignment is represented
                 as a tuple of aligned strings
        """
        
        
        return self.align()
        """
        return [
            ('ADMI-NS', 'ADMIRES'), ('ADMIN-S', 'ADMIRES')
        ]
        """

    def get_score_matrix(self):
        """
        :return: matrix built during the alignment process as a list of lists
        """
        """
        return [
            [0, -1, -2, -3, -4, -5, -6],
            [-1, 1, 0, -1, -2, -3, -4],
            [-2, 0, 2, 1, 0, -1, -2],
            [-3, -1, 1, 3, 2, 1, 0],
            [-4, -2, 0, 2, 4, 3, 2],
            [-5, -3, -1, 1, 3, 4, 3],
            [-6, -4, -2, 0, 2, 3, 4],
            [-7, -5, -3, -1, 1, 2, 4]
        ]
        """
        return self.score_matrix

"""
input = {
        "strings": ["CYVPST", "WYVPST"],
        "gap_penalty": -6,
        "matrix": "identity",
        "best_score": 0,
        "number_of_alignments": 2,
        "alignments": [["AVNCCEGQHI", "ARN-DE-Q--"], ["AVNCCEGQHI", "ARND-E-Q--"]]
    }
g =  GlobalAlignment(input["strings"][0], input["strings"][1], input["gap_penalty"], MATRICES[input['matrix']])
"""
