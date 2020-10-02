import numpy as np


"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix or array provided as
           parameters. Further, use those indices when generating or returning
           any matrices or arrays. Failure to do so will most likely result in
           not passing the tests.
EXAMPLE: To access the substitution frequency from alanine 'A' to proline 'P'
         in the bg_matrix use bg_matrix[AA_TO_INT['A'], AA_TO_INT['P']].
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}
GAP_INDEX = AA_TO_INT['-']


class MSA:

    def __init__(self, sequences):
        """
        Initialize the MSA class with the provided list of sequences. Check the
        sequences for correctness. Pre-calculate any statistics you seem fit.

        :param sequences: List containing the MSA sequences.
        """
        if not sequences:
            raise TypeError
        ref_len = len(sequences[0])
        for seq in sequences:
            if len(seq) != ref_len:
                raise TypeError
            valid_char = set([ch in ALPHABET for ch in seq])
            if len(valid_char) == 2:
                raise TypeError
        self.sequences = sequences
        self.seq_mat = np.array(self.sequences)

    def get_pssm(self, *, bg_matrix=None, beta=10, use_sequence_weights=False,
                 redistribute_gaps=False, add_pseudocounts=False):
        """
        Return a PSSM for the underlying MSA. Use the appropriate refinements 
        according to the parameters. If no bg_matrix is specified, use uniform 
        background and pair frequencies.
        Every row in the resulting PSSM corresponds to a non-gap position in 
        the primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).

        :param bg_matrix: Amino acid pair frequencies as numpy array (20, 20).
                          Access the matrix using the indices from AA_TO_INT.
        :param beta: Beta value (float) used to weight the pseudocounts 
                     against the observed amino acids in the MSA.
        :param use_sequence_weights: Calculate and apply sequence weights.
        :param redistribute_gaps: Redistribute the gaps according to the 
                                  background frequencies.
        :param add_pseudocounts: Calculate and add pseudocounts according 
                                 to the background frequencies.

        :return: PSSM as numpy array of shape (L x 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        BG_FREQ = 0.05
        size = self.get_size()
        seqs = self.sequences
        pssm = np.zeros((size[1], 21))
        if not use_sequence_weights:
            for col in range(size[1]):
                for i in range(len(seqs)):
                    if not redistribute_gaps and seqs[i][col] != '-':
                        indx = AA_TO_INT[seqs[i][col]]
                        pssm[col][indx] += 1
                    elif redistribute_gaps:
                        indx = AA_TO_INT[seqs[i][col]]
                        pssm[col][indx] += 1
        else:
            seq_weights = self.get_sequence_weights()
            for col in range(size[1]):
                for i in range(len(seqs)):
                    
                    indx = AA_TO_INT[seqs[i][col]]
                    pssm[col][indx] += seq_weights[i]
            if not redistribute_gaps:
                pssm = np.delete(pssm, -1, axis=1)  
        if redistribute_gaps:
            for i in range(pssm.shape[0]):
                if bg_matrix is None and pssm[i][-1] > 0:
                    mul_factor = BG_FREQ * pssm[i][-1]
                    pssm[i][:] = pssm[i][:] + mul_factor 
                elif bg_matrix is not None and pssm[i][-1] > 0:
                    
                    for j in range(pssm.shape[1] - 1):
                        bg = np.sum(bg_matrix[AA_TO_INT[ALPHABET[j]]])
                        mul_factor = bg * pssm[i][-1]
                        pssm[i][j] = pssm[i][j] + mul_factor
                    """
                    mul_factor = BG_FREQ * pssm[i][-1]
                    pssm[i][:] = pssm[i][:] + mul_factor
                    """
            pssm = np.delete(pssm, -1, axis=1)   
        # Pseudocount
        if add_pseudocounts:
            ps_mat = np.zeros((size[1], 20))  
            if bg_matrix is not None:   
                for i in  range(pssm.shape[0]):
                
                    for j in range(pssm.shape[1]):
                        ps_count = 0
                        for k in range(pssm.shape[1]):
                            
                            g = bg_matrix[AA_TO_INT[ALPHABET[j]]][AA_TO_INT[ALPHABET[k]]]
                            g_calc = g * pssm[i][k] / np.sum(bg_matrix[AA_TO_INT[ALPHABET[k]]])
                            ps_count += g_calc
                        ps_mat[i, j] = ps_count
            else:
                #g = bg_matrix[AA_TO_INT[ALPHABET[0]]]
                """
                pssm = np.delete(pssm, -1, axis=1) 
                ps_mat[:,:] = (pssm[:, :] + 1)* BG_FREQ
                
                
                """
                pssm = np.delete(pssm, -1, axis=1) 
                for i in  range(pssm.shape[0]):
                    for j in range(pssm.shape[1]):
                        ps_count = 0
                        for k in range(pssm.shape[1]):
                                
                            #g = bg_matrix[AA_TO_INT[ALPHABET[j]]][AA_TO_INT[ALPHABET[k]]]
                            #g_calc = pssm[i][k] * pssm[i][k] / BG_FREQ
                            g_calc = pssm[i][k] * BG_FREQ
                            ps_count += g_calc
                        #ps_mat[i, j] = pssm[i,j] + ps_count
                        ps_mat[i, j] = ps_count
                #pssm = ps_mat
                
                
                """
                max_indx = np.where(pssm[i][:] = max(pssm[i][:]))
                g_self = bg_matrix[AA_TO_INT[ALPHABET[max_indx]]][AA_TO_INT[ALPHABET[max_indx]]]
                calc = g_self * pssm[i][indx] / BG_FREQ
                ps_mat[i][max_indx] = calc
                """
            if True:
                N = self.get_number_of_observations()
                alpha = N - 1
                temp_pssm = np.zeros((size[1], 20))   
                for i in  range(pssm.shape[0]):
                    for j in range(pssm.shape[1]):
                        adj_f = ((alpha * pssm[i][j]) + (beta * ps_mat[i][j]))/(alpha + beta)
                        temp_pssm[i][j] = adj_f
                pssm = temp_pssm

        # Normalize
        for i in range(pssm.shape[0]):
            sum = np.sum(pssm[i, :])
            pssm[i, :] = pssm[i, :]/sum
        if bg_matrix is None:
            pssm[:, :] = pssm[:, :]/BG_FREQ
        else:
            for i in range(pssm.shape[0]):
                if redistribute_gaps or use_sequence_weights or add_pseudocounts:
                    range_num = pssm.shape[1]
                else:
                    range_num = pssm.shape[1] - 1
                for j in range(range_num):
                    bg = np.sum(bg_matrix[AA_TO_INT[ALPHABET[j]]])
                    pssm[i, j] = pssm[i, j]/bg
        pssm[:, :] = 2*np.log2(pssm[:,:])
        pssm[np.where(np.isneginf(pssm))] = -20
        #Every row in the resulting PSSM corresponds to a non-gap position in 
        #the primary sequence of the MSA (i.e. the first one).
        pri_seq = self.sequences[0]
        indxs = [i for i, letter in enumerate(pri_seq) if letter == '-']
        pssm = np.delete(pssm, indxs, axis=0)
        if not redistribute_gaps and not use_sequence_weights and not add_pseudocounts:
            pssm = np.delete(pssm, -1, axis=1)
        return np.rint(pssm).astype(np.int64)

    def get_size(self):
        """
        Return the number of sequences in the MSA and the MSA length, i.e.
        the number of columns in the MSA. This includes gaps.

        :return: Tuple of two integers. First element is the number of
                 sequences in the MSA, second element is the MSA length.
        """
        return (len(self.sequences), len(self.sequences[0]))

    def get_primary_sequence(self):
        """
        Return the primary sequence of the MSA. In this exercise, the primary
        sequence is always the first sequence of the MSA. The returned 
        sequence must NOT include gap characters.

        :return: String containing the ungapped primary sequence.
        """
        pri_seq = self.sequences[0].replace("-", "")
        return pri_seq

    def get_sequence_weights(self):
        """
        Return the calculated sequence weights for all sequences in the MSA.
        The order of weights in the array must be equal to the order of the
        sequences in the MSA.

        :return: Numpy array (dtype=numpy.float64) containing the weights for
                 all sequences in the MSA.
        """
        #weights = np.zeros(10)
        size = self.get_size()
        seqs = self.sequences
        # plus 1 for storing r
        w_mat = np.zeros((size[1], size[0] + 1))
        for i in range(size[1]):
            aaciddict = {}
            # Create the amino acid list at a particular column
            for j in range(size[0]):
                if seqs[j][i] not in aaciddict:
                    aaciddict[seqs[j][i]] = 1
                else:
                    aaciddict[seqs[j][i]] += 1
            # Find the number of unique aa
            r = len(set(aaciddict.keys()))
            w_mat[i][-1] = r
            for a in range(size[0]):
                freq = aaciddict[seqs[a][i]]
                w_mat[i][a] = 1/(r*freq) 
        tot_sum = []
        for i in range(size[0]):
            sum = 0
            for j in range(size[1]):
                
                if w_mat[j][-1] > 1:
                    sum += w_mat[j][i]
            tot_sum.append(sum)
        self.weight_mat = w_mat
        return np.array(tot_sum).astype(np.float64)

    def get_number_of_observations(self):
        """
        Return the estimated number of independent observations in the MSA.

        :return: Estimate of independent observation (dtype=numpy.float64).
        """
        num_obs = np.int64(-1)
        w_mat = self.weight_mat
        sum_r = np.sum(w_mat[:, -1], axis=0)
        num_obs = sum_r/(np.shape(w_mat)[0])
        return num_obs.astype(np.float64)
"""
from pathlib import Path
import json

test_json = 'pssm_test.json'
relative_path = "/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise4-ge73tag/tests/"

with Path(relative_path, test_json).open('r') as json_file:
    json_data = json.load(json_file)
msa = json_data['msa_sequences']
inv_msa = json_data['invalid_msa']
bg_matrix = json_data['bg_matrix']
msa = MSA(msa)
msa.get_sequence_weights()
msa.get_number_of_observations()

#msa.get_pssm(bg_matrix=bg_matrix, redistribute_gaps=True)
msa.get_pssm(add_pseudocounts=True)
"""