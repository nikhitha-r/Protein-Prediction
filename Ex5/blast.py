import numpy as np

from pathlib import Path
import re
from time import time as ttime
import json
"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix (PSSM or substitution
           matrix) parameters. Failure to do so will most likely result in not
           passing the tests.
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}


class BlastDB:

    def __init__(self):
        """
        Initialize the BlastDB class.
        """
        self.blastdb = []
        self.match_seqs = {}
        

    def add_sequence(self, sequence):
        """
        Add a sequence to the database.

        :param sequence: a protein sequence (string).
        """
        self.blastdb.append(sequence)

    def get_sequences(self, word):
        """
        Return all sequences in the database containing a given word.

        :param word: a word (string).

        :return: List with sequences.
        """
        match_seqs = []
        for seq in self.blastdb:
            if word in seq:
                match_seqs.append(seq)
        return match_seqs

    def get_seq_two_hit(self, word, word_pos_query):
        for seq in self.blastdb:
            if word in seq:
                stri = "(?=" + word + ")"
                strt_indxs = [m.start() for m in re.finditer(stri, seq)]
                if seq not in self.match_seqs:
                    
                    self.match_seqs[seq] = [[word_pos_query[word], strt_indxs, word]]
                else:
                    self.match_seqs[seq].append([word_pos_query[word], strt_indxs, word])
            
            #self.match_seqs = {k: v for k, v in sorted(self.match_seqs.items(), key=lambda item: item[0][0])}
            """
            if seq in self.match_seqs and len(self.match_seqs[seq]) > 1:
                self.match_seqs[seq].sort(key = lambda test_list: test_list[0][0])
            """
        return self.match_seqs

    def get_db_stats(self):
        """
        Return some database statistics:
            - Number of sequences in database
            - Number of different words in database
            - Average number of words per sequence (rounded to nearest int)
            - Average number of sequences per word (rounded to nearest int)

        :return: Tuple with four integer numbers corrsponding to the mentioned
                 statistics (in order of listing above).
        """
        num_seq = len(self.blastdb)
        word_list = []
        len_avg_word_per_seq_list = []
        word_dict = {}
        for seq in self.blastdb:
           # for i in range(3):
            """
            if True:
                start = 0
                for j in range(start, len(seq), 3):
                    word = seq[j:j+3]
                    if word not in word_list and len(word) == 3:
                        word_list.append(word)
            """
            split_seqs = []
            for i in range(3):
                split_seq = [seq[start:start+3] for start in range(i, len(seq), 3)]
                if len(split_seq[-1]) != 3:
                    del split_seq[-1]
                
                #split_seq = list(set(split_seq))
                split_seqs = split_seqs + split_seq
                word_list  = word_list + split_seq
                word_list = list(set(word_list))
            split_seqs = set(split_seqs)
            len_avg_word_per_seq_list.append(len(split_seqs))
            split_seqs = set(split_seqs)
            for word in split_seqs:
                if word not in word_dict:
                    word_dict[word] = 1
                else:
                    word_dict[word] += 1
        num_diff_words = len(word_list)
        avg_words = np.round(np.mean(len_avg_word_per_seq_list)).astype(np.int)
        avg_seq = np.round(np.mean(list(word_dict.values()))).astype(np.int)
        return (num_seq, num_diff_words, avg_words, avg_seq)


class Blast:

    def __init__(self, substitution_matrix):
        """
        Initialize the Blast class with the given substitution_matrix.

        :param substitution_matrix: 20x20 amino acid substitution score matrix.
        """
        self.sub_matrix = substitution_matrix
        

    def get_words(self, *, sequence=None, pssm=None, T=11):
        """
        Return all words with score >= T for given protein sequence or PSSM.
        Only a sequence or PSSM will be provided, not both at the same time.
        A word may only appear once in the list.

        :param sequence: a protein sequence (string).
        :param pssm: a PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.

        :return: List of unique words.
        """
        """
        split_seqs = []
        for i in range(3):
            split_seq = [sequence[start:start+3] for start in range(i, len(sequence), 3)]
            if len(split_seq[-1]) != 3:
                del split_seq[-1]
            
            split_seq = list(set(split_seq))
            split_seqs = split_seqs + split_seq
        
        uniq_words = []
        for word in split_seqs:
            score = 0
            for ch in word:
                score += AA_TO_INT[ch]
            if score >= T:
                uniq_words.append(word)
        print(split_seqs)
        print(uniq_words)
        """
        uniq_words = []
        self.word_pos = {}
        if pssm is None:
            for i in range(len(sequence)):
                word = sequence[i:i+3]
                if len(word) == 3:
                    for a in ALPHABET:
                        for b in ALPHABET:
                            for c in ALPHABET:
                                score = self.sub_matrix[AA_TO_INT[a]][AA_TO_INT[sequence[i]]] \
                                        + self.sub_matrix[AA_TO_INT[b]][AA_TO_INT[sequence[i + 1]]] \
                                        + self.sub_matrix[AA_TO_INT[c]][AA_TO_INT[sequence[i + 2]]]
                                if score >= T and (a+b+c) not in uniq_words:
                                    uniq_words.append(a+b+c)
                                    self.word_pos[a+b+c] = [i]
                                elif score >= T and (a+b+c) in uniq_words:
                                    self.word_pos[a+b+c].append(i)
                                #if score >= T and word not in uniq_words:
                                #    uniq_words.append(word)
                                

        else:
            print(len(pssm))
            for i in range(len(pssm) - 2):
                for a in range(20):
                    for b in range(20):
                        for c in range(20):
                            #score = pssm[AA_TO_INT[a]][0] + pssm[AA_TO_INT[b]][1] + pssm[AA_TO_INT[c]][2]
                            score = pssm[i][a] + pssm[i+1][b] + pssm[i+2][c]
                            #print(score)
                            word = INT_TO_AA[a] + INT_TO_AA[b] + INT_TO_AA[c]
                            if score >= T and word not in uniq_words:
                                uniq_words.append(word)
                                self.word_pos[word] = [i]
                            elif score >= T and word in uniq_words:
                                self.word_pos[word].append(i)
        #print(uniq_words)
        self.uniq_words = uniq_words            
        return uniq_words

    def search_one_hit(self, blast_db, *, query=None, pssm=None, T=13, X=5, S=30):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.

        :return: dictionary of target sequences and list of HSP tuples.
        """
        #print("ADADADADAADAADA", query)
        if pssm is None:
            words = self.get_words(sequence=query, pssm=None, T=T)
        else:
            words = self.get_words(sequence=query, pssm=pssm, T=T) 
            self.sub_matrix = pssm
        if pssm is None:
            result = {}
            res_seqs = {}
            for word in words:
                #######################

                #word = 'ERP'
                ########################
                seqs = blast_db.get_sequences(word)
                #print("##############", seqs)
                #strt_indxs_qu = [m.start() for m in re.finditer(word, query)][0]
                strt_indxs_qus = self.word_pos[word]
                
                for seq in seqs:
                    
                    ##############
                    if True:
                    #if seq in list(blast_ex.keys()) and (seq == 'MLTMSVTLSPLRSQGPDPMATDASPMAINMTPTVEQEEGEGEEAVKAIDAEQQYGKPPPLHTAADWKIVLHLPEIETWLRMTSERVRDLTYSVQQDADSKHVDVHLVQLKDICEDISDHVEQIHALLETEFSLKLLSYSVNVIVDIHAVQLLWHQLRVSVLVLRERILQGLQDANGNYTRQTDILQAFSEETTEGRLDSLTEVDDSGQLTIKCSQDYLSLDCGITAFELSDYSPSEDLLGGLGDMTTSQAKTKSFDSWSYSEMEKEFPELIRSVGLLTVATEPVPSSCGEANEDSSQASLSDDHKGEHGEDGAPVPGQQLDSTVGMSSLDGTLANAAEHPSETAKQDSTSSPQLGAKKTQPGPCEITTPKRSIRDCFNYNEDSPTQPTLPKRGLFLKETQKNERKGSDRKGQVVDLKPELSRSTPSLVDPPDRSKLCLVLQSSYPSSPSAASQSYECLHKVGLGNLENIVRSHIKEISSSLGRLTDCHKEKLRLKKPHKTLAEVSLCRIPKQGGGSGKRSESTGSSAGPSMVSPGAPKATMRPETDSASTASGGLCHQRNRSGQLPVQSKASSSPPCSHSSESSLGSDSIKSPVPLLSKNKSQKSSPPAPCHATQNGQVVEAWYGSDEYLALPSHLKQTEVLALKLESLTKLLPQKPRGETIQDIDDWELSEMNSDSEIYPTYHIKKKHTRLGTVSPSSSSDIASSLGESIESGPLSDILSDEDLCLPLSSVKKFTDEKSERPSSSEKNESHSATRSALIQKLMHDIQHQENYEAIWERIEGFVNKLDEFIQWLNEAMETTENWTPPKAETDSLRLYLETHLSFKLNVDSHCALKEAVEEEGHQLLELVVSHKAGLKDTLRMIASQWKELQRQIKRQHSWILRALDTIKAEILATDVSVEDEEGTGSPKAEVQLCHLETQRDAVEQMSLKLYSEQYTSGSKRKEEFANMSKAHAEGSNGLLDFDSEYQELWDWLIDMESLVMDSHDLMMSEEQQQHLYKRYSVEMSIRHLKKSELLSKVEALKKGGLSLPDDILEKVDSINEKWELLGKTLREKIQDTIAGHSGSGPRDLLSPESGSLVRQLEVRIKELKRWLRDTELFIFNSCLRQEKEGTSAEKQLQYFKSLCREIKQRRRGVASILRLCQHLLDDRDTCNLNADHQPMQLIIVNLERRWEAIVMQAVQWQTRLQKKMGKESETLNVIDPGLMDLNGMSEDALEWDETDISNKLISVHEESNDLDQDPEPMLPAVKLEETHHKDSGYEEEAGDCGGSPYTSNITAPSSPHIYQVYSLHNVELHEDSHTPFLKSSPKFTGTTQPTVLTKSLSKDSSFSSTKSLPDLLGGSGLVRPYSCHSGDLSQNSGSESGIVSEGDNEMPTNSDMSLFSMVDGSPSNPETEHPDPQMGDAANVLEQKFKDNGESIKLSSVSRASVSPVGCVNGKAGDLNSVTKHTADCLGEELQGKHDVFTFYDYSYLQGSKLKLPMIMKQPQSEKAHVEDPLLGGFYFDKKSCKAKHQASESQPDAPPHERILASAPHEMGRSAYKSSDIEKTFTGIQSARQLSLLSRSSSVESLSPGGDLFGLGIFKNGSDSLQRSTSLESWLTSYKSNEDLFSCHSSGDISVSSGSVGELSKRTLDLLNRLENIQSPSEQKIKRSVSDMTLQSSSQKMPFAGQMSLDVASSINEDSPASLTELSSSDELSLCSEDIVLHKNKIPESNASFRKRLNRSVADESDVNVSMIVNVSCTSACTDDEDDSDLLSSSTLTLTEEELCLKDEDDDSSIATDDEIYEESNLMSGLDYIKNELQTWIRPKLSLTREKKRSGVTDEIKVNKDGGGNEKANPSDTLDIEALLNGSIRCLSENNGNGKTPPRTHGSGTKGENKKSTYDVSKDPHVADMENGNIESTPEREREKPQGLPEVSENLASNVKTISESELSEYEAVMDGSEDSSVARKEFCPPNDRHPPQMGPKLQHPENQSGDCKPVQNPCPGLLSEAGVGSRQDSNGLKSLPNDAPSGARKPAGCCLLEQNETEESASISSNASCCNCKPDVFHQKDDEDCSVHDFVKEIIDMASTALKSKSQPESEVAAPTSLTQIKEKVLEHSHRPIHLRKGDFYSYLSLSSHDSDCGEVTNYIDEKSSTPLPPDAVDSGLDDKEDMDCFFEACVEDEPVNEEAGLPGALPNESAIEDGAEQKSEQKTASSPVLSDKTDLVPLSGLSPQKGADDAKEGDDVSHTSQGCAESTEPTTPSGKANAEGRSRMQGVSATPEENAASAKPKIQAFSLNAKQPKGKVAMRYPSPQTLTCKEKLVNFHEDRHSNMHR'):
                    #############
                        stri = "(?=" + word + ")"
                        strt_indxs = [m.start() for m in re.finditer(stri, seq)]
                        print(len(result))

                        #strt_indx = seq.index(word)
                        # Calculate the initial score of the word
                        for strt_indx in strt_indxs:
                            for strt_indxs_qu in strt_indxs_qus:
                                temp_seq = word
                                diff_indx = strt_indxs_qu - strt_indx
                                score_list = []
                                seq_score = {}
                                print(AA_TO_INT[query[strt_indx + diff_indx]])
                                score = self.sub_matrix[AA_TO_INT[query[strt_indx + diff_indx]]][AA_TO_INT[seq[strt_indx]]] + \
                                        self.sub_matrix[AA_TO_INT[query[strt_indx + 1 + diff_indx]]][AA_TO_INT[seq[strt_indx + 1]]] + \
                                        self.sub_matrix[AA_TO_INT[query[strt_indx + 2 + diff_indx]]][AA_TO_INT[seq[strt_indx + 2]]]
                                print(score)
                                score_list.append(score)
                                seq_score[temp_seq] = [score, strt_indx]

                                

                                if (score - (max(score_list))) < X:
                                    backtrack = False
                                    for i in range(strt_indx + 3, len(seq)):
                                        if (i + diff_indx) < len(query):
                                            score += self.sub_matrix[AA_TO_INT[query[i + diff_indx]]][AA_TO_INT[seq[i]]]                        
                                            #if (score - (max(score_list))) < X and (i + diff_indx) != (len(query) - 1):
                                            #if ((max(score_list)) - score) < X and (i + diff_indx) != (len(query) - 1) and score < S:
                                            if ((max(score_list)) - score) < X and (i + diff_indx) != (len(query) - 1):
                                                temp_seq += seq[i]
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, False]
                                            elif (i + diff_indx) == (len(query) - 1) or (max(score_list) - score) >= X:
                                                temp_seq += seq[i]
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, False]
                                                backtrack = True
                                                
                                                """
                                                sc = [item[0] for item in list(seq_score.values())]
                                                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                                max_value_seq = min(max_value_seqs, key=len)
                                                max_strt_idx = seq_score[max_value_seq][1]
                                                temp_seq = max_value_seq
                                                score = max(sc)
                                                """
                                                break
                                            """
                                            elif score >= S:
                                                temp_seq += seq[i]
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, False]
                                                break
                                            """
                                        else:
                                            backtrack = True
                                            break

                                    if backtrack or i == (len(seq) - 1) or (strt_indx + 3 == len(seq)):
                                        sc = [item[0] for item in list(seq_score.values())]
                                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                        max_value_seq = min(max_value_seqs, key=len)
                                        max_strt_idx = seq_score[max_value_seq][1]
                                        temp_seq = max_value_seq
                                        score = max(sc)
                                        #for i in range(strt_indx - 1, 0 , -1):
                                        for i in range(strt_indx - 1, -1 , -1):
                                            score += self.sub_matrix[AA_TO_INT[query[i + diff_indx]]] [AA_TO_INT[seq[i]]]                        
                                            #if ((max(score_list)) - score) < X and (i + diff_indx) != 0 and score < S:
                                            if ((max(score_list)) - score) < X and (i + diff_indx) != 0:
                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, True]
                                                strt_b_indx = i
                                            elif (i + diff_indx) == 0 or (max(score_list) - score) >= X:

                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, True]

                                                backtrack = True
                                                """
                                                sc = [item[0] for item in list(seq_score.values())]
                                                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                                max_value_seq = min(max_value_seqs, key=len)
                                                max_strt_idx = seq_score[max_value_seq][1]
                                                temp_seq = max_value_seq
                                                score = max(sc)
                                                """
                                                """
                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i]
                                                
                                                
                                            
                                                sc = [item[0] for item in list(seq_score.values())]
                                                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                                max_value_seq = max(max_value_seqs, key=len)
                                                max_strt_idx = seq_score[max_value_seq][1]
                                                """
                                                break
                                            elif score >= S:
                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i]
                                                break
                                    if backtrack or (backtrack and i == 0) or (not backtrack and (i == len(seq) - 1)):
                                        sc = [item[0] for item in list(seq_score.values())]
                                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                        max_value_seq = min(max_value_seqs, key=len)
                                        max_strt_idx = seq_score[max_value_seq][1]
                                        temp_seq = max_value_seq
                                        score = max(sc) 
                                    print(len(result))

                                    if seq_score[temp_seq][0] < S:
                                        continue
                                    else:
                                        
                                        sc = [item[0] for item in list(seq_score.values())]
                                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                        max_value_seq = min(max_value_seqs, key=len)
                                        temp_seq = max_value_seq
                                        #tar_idx = min(strt_indx, strt_b_indx)
                                        #res = [min(strt_indx, strt_b_indx) + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0])]
                                        if seq_score[max_value_seq][2]:
                                            tar_idx = seq_score[max_value_seq][1]
                                        else:
                                            tar_idx = seq_score[max_value_seq][1] - len(max_value_seq) + 1
                                        if (tar_idx + diff_indx) >= 0:
                                            res = (tar_idx + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0]))
                                            target = seq[:tar_idx] + max_value_seq + seq[tar_idx+len(max_value_seq):]
                                            #target[tar_idx:(tar_idx+len(max_value_seq))] = max_value_seq
                                            if target not in result:
                                                result[target] = [res]
                                                mod_res = [res, temp_seq]
                                                res_seqs[target] = [mod_res]
                                            else:
                                                present = False
                                                res_seq = result[target]
                                                for res_list in res_seq:
                                                    indx = res_list
                                                    #if (res[1] + res[2] == indx[1] + indx[2]) or (res[0]) in range(indx[0], indx[0] + indx[2]) or (indx[0]) in range(res[0], res[0] + res[2]):
                                                    if (res[1] + res[2] == indx[1] + indx[2]) and (res[0] + res[2] == indx[0] + indx[2]):
                                                        present = True
                                                        if res[3] >= indx[3]:
                                                            if res[0] < indx[0] and indx in result[target]:
                                                                print(indx)
                                                                result[target].remove(indx)
                                                        # if (res[3] != indx[3] and (res[0] + res[2] not in [31] or res[3] == 43)) or res[0] + res[2] in [32, 64, 21, 23] or  res[1] + res[2] in [47]:
                                                            #result[target].remove(indx)
                                                            if res not in result[target]:
                                                                result[target].append(res)
                                                if not present:
                                                    old_list = result[target]
                                                    if res not in old_list:
                                                        result[target].append(res)
                                                        mod_res = [res, temp_seq]
                                                        res_seqs[target].append(mod_res)
                                    #result[seq] = [[min(strt_indx, strt_b_indx) + diff_indx, min(strt_indx, strt_b_indx),  len(temp_seq), int(seq_score[temp_seq][0]])]

        
        #PSSM#########
        if pssm is not None:
            words = self.get_words(sequence=query, pssm=pssm, T=T)
            
            result = {}
            res_seqs = {}
            print("@@@@@@@@@@@@")
            for word in words:
                #######################

                #word = 'ERP'
                ########################
                seqs = blast_db.get_sequences(word)
                #print("##############", seqs)
                #strt_indxs_qu = [m.start() for m in re.finditer(word, query)][0]
                strt_indxs_qus = self.word_pos[word]
                
                for seq in seqs:
                    
                    ##############
                    if True:
                    #if seq in list(blast_ex.keys()) and (seq == 'MLTMSVTLSPLRSQGPDPMATDASPMAINMTPTVEQEEGEGEEAVKAIDAEQQYGKPPPLHTAADWKIVLHLPEIETWLRMTSERVRDLTYSVQQDADSKHVDVHLVQLKDICEDISDHVEQIHALLETEFSLKLLSYSVNVIVDIHAVQLLWHQLRVSVLVLRERILQGLQDANGNYTRQTDILQAFSEETTEGRLDSLTEVDDSGQLTIKCSQDYLSLDCGITAFELSDYSPSEDLLGGLGDMTTSQAKTKSFDSWSYSEMEKEFPELIRSVGLLTVATEPVPSSCGEANEDSSQASLSDDHKGEHGEDGAPVPGQQLDSTVGMSSLDGTLANAAEHPSETAKQDSTSSPQLGAKKTQPGPCEITTPKRSIRDCFNYNEDSPTQPTLPKRGLFLKETQKNERKGSDRKGQVVDLKPELSRSTPSLVDPPDRSKLCLVLQSSYPSSPSAASQSYECLHKVGLGNLENIVRSHIKEISSSLGRLTDCHKEKLRLKKPHKTLAEVSLCRIPKQGGGSGKRSESTGSSAGPSMVSPGAPKATMRPETDSASTASGGLCHQRNRSGQLPVQSKASSSPPCSHSSESSLGSDSIKSPVPLLSKNKSQKSSPPAPCHATQNGQVVEAWYGSDEYLALPSHLKQTEVLALKLESLTKLLPQKPRGETIQDIDDWELSEMNSDSEIYPTYHIKKKHTRLGTVSPSSSSDIASSLGESIESGPLSDILSDEDLCLPLSSVKKFTDEKSERPSSSEKNESHSATRSALIQKLMHDIQHQENYEAIWERIEGFVNKLDEFIQWLNEAMETTENWTPPKAETDSLRLYLETHLSFKLNVDSHCALKEAVEEEGHQLLELVVSHKAGLKDTLRMIASQWKELQRQIKRQHSWILRALDTIKAEILATDVSVEDEEGTGSPKAEVQLCHLETQRDAVEQMSLKLYSEQYTSGSKRKEEFANMSKAHAEGSNGLLDFDSEYQELWDWLIDMESLVMDSHDLMMSEEQQQHLYKRYSVEMSIRHLKKSELLSKVEALKKGGLSLPDDILEKVDSINEKWELLGKTLREKIQDTIAGHSGSGPRDLLSPESGSLVRQLEVRIKELKRWLRDTELFIFNSCLRQEKEGTSAEKQLQYFKSLCREIKQRRRGVASILRLCQHLLDDRDTCNLNADHQPMQLIIVNLERRWEAIVMQAVQWQTRLQKKMGKESETLNVIDPGLMDLNGMSEDALEWDETDISNKLISVHEESNDLDQDPEPMLPAVKLEETHHKDSGYEEEAGDCGGSPYTSNITAPSSPHIYQVYSLHNVELHEDSHTPFLKSSPKFTGTTQPTVLTKSLSKDSSFSSTKSLPDLLGGSGLVRPYSCHSGDLSQNSGSESGIVSEGDNEMPTNSDMSLFSMVDGSPSNPETEHPDPQMGDAANVLEQKFKDNGESIKLSSVSRASVSPVGCVNGKAGDLNSVTKHTADCLGEELQGKHDVFTFYDYSYLQGSKLKLPMIMKQPQSEKAHVEDPLLGGFYFDKKSCKAKHQASESQPDAPPHERILASAPHEMGRSAYKSSDIEKTFTGIQSARQLSLLSRSSSVESLSPGGDLFGLGIFKNGSDSLQRSTSLESWLTSYKSNEDLFSCHSSGDISVSSGSVGELSKRTLDLLNRLENIQSPSEQKIKRSVSDMTLQSSSQKMPFAGQMSLDVASSINEDSPASLTELSSSDELSLCSEDIVLHKNKIPESNASFRKRLNRSVADESDVNVSMIVNVSCTSACTDDEDDSDLLSSSTLTLTEEELCLKDEDDDSSIATDDEIYEESNLMSGLDYIKNELQTWIRPKLSLTREKKRSGVTDEIKVNKDGGGNEKANPSDTLDIEALLNGSIRCLSENNGNGKTPPRTHGSGTKGENKKSTYDVSKDPHVADMENGNIESTPEREREKPQGLPEVSENLASNVKTISESELSEYEAVMDGSEDSSVARKEFCPPNDRHPPQMGPKLQHPENQSGDCKPVQNPCPGLLSEAGVGSRQDSNGLKSLPNDAPSGARKPAGCCLLEQNETEESASISSNASCCNCKPDVFHQKDDEDCSVHDFVKEIIDMASTALKSKSQPESEVAAPTSLTQIKEKVLEHSHRPIHLRKGDFYSYLSLSSHDSDCGEVTNYIDEKSSTPLPPDAVDSGLDDKEDMDCFFEACVEDEPVNEEAGLPGALPNESAIEDGAEQKSEQKTASSPVLSDKTDLVPLSGLSPQKGADDAKEGDDVSHTSQGCAESTEPTTPSGKANAEGRSRMQGVSATPEENAASAKPKIQAFSLNAKQPKGKVAMRYPSPQTLTCKEKLVNFHEDRHSNMHR'):
                    #############
                        stri = "(?=" + word + ")"
                        strt_indxs = [m.start() for m in re.finditer(stri, seq)]
                        
                        #strt_indx = seq.index(word)
                        # Calculate the initial score of the word
                        for strt_indx in strt_indxs:
                            for strt_indxs_qu in strt_indxs_qus:
                                temp_seq = word
                                diff_indx = strt_indxs_qu - strt_indx
                                score_list = []
                                seq_score = {}
                                """
                                score = self.sub_matrix[AA_TO_INT[query[strt_indx + diff_indx]]][AA_TO_INT[seq[strt_indx]]] + \
                                        self.sub_matrix[AA_TO_INT[query[strt_indx + 1 + diff_indx]]][AA_TO_INT[seq[strt_indx + 1]]] + \
                                        self.sub_matrix[AA_TO_INT[query[strt_indx + 2 + diff_indx]]][AA_TO_INT[seq[strt_indx + 2]]]
                                """
                                print("@@@@@@@@@@@@")
                                score = pssm[strt_indxs_qu][AA_TO_INT[seq[strt_indx]]] + \
                                        pssm[strt_indxs_qu + 1][AA_TO_INT[seq[strt_indx + 1]]] + \
                                        pssm[strt_indxs_qu + 2][AA_TO_INT[seq[strt_indx + 2]]] 
                                score_list.append(score)
                                print(score)
                                seq_score[temp_seq] = [score, strt_indx]
                                if (score - (max(score_list))) < X:
                                    backtrack = False
                                    for i in range(strt_indx + 3, len(seq)):
                                        if (i + diff_indx) < len(pssm):
                                            score += self.sub_matrix[i + diff_indx][AA_TO_INT[seq[i]]] 
                                                                   
                                            #if (score - (max(score_list))) < X and (i + diff_indx) != (len(query) - 1):
                                            #if ((max(score_list)) - score) < X and (i + diff_indx) != (len(query) - 1) and score < S:
                                            if ((max(score_list)) - score) < X and (i + diff_indx) != (len(pssm) - 1):
                                                temp_seq += seq[i]
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, False]
                                            elif (i + diff_indx) == (len(pssm) - 1) or (max(score_list) - score) >= X:
                                                temp_seq += seq[i]
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, False]
                                                backtrack = True
                                                
                                                """
                                                sc = [item[0] for item in list(seq_score.values())]
                                                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                                max_value_seq = min(max_value_seqs, key=len)
                                                max_strt_idx = seq_score[max_value_seq][1]
                                                temp_seq = max_value_seq
                                                score = max(sc)
                                                """
                                                break
                                            """
                                            elif score >= S:
                                                temp_seq += seq[i]
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, False]
                                                break
                                            """
                                        else:
                                            backtrack = True
                                            break
                                    print("@@@@@@@@")
                                    if backtrack or i == (len(seq) - 1) or (strt_indx + 3 == len(seq)):
                                        sc = [item[0] for item in list(seq_score.values())]
                                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                        max_value_seq = min(max_value_seqs, key=len)
                                        max_strt_idx = seq_score[max_value_seq][1]
                                        temp_seq = max_value_seq
                                        score = max(sc)
                                        #for i in range(strt_indx - 1, 0 , -1):
                                        for i in range(strt_indx - 1, -1 , -1):
                                            score += self.sub_matrix[i + diff_indx] [AA_TO_INT[seq[i]]]                        
                                            #if ((max(score_list)) - score) < X and (i + diff_indx) != 0 and score < S:
                                            if ((max(score_list)) - score) < X and (i + diff_indx) != 0:
                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, True]
                                                strt_b_indx = i
                                            elif (i + diff_indx) == 0 or (max(score_list) - score) >= X:

                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i, True]

                                                backtrack = True
                                                """
                                                sc = [item[0] for item in list(seq_score.values())]
                                                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                                max_value_seq = min(max_value_seqs, key=len)
                                                max_strt_idx = seq_score[max_value_seq][1]
                                                temp_seq = max_value_seq
                                                score = max(sc)
                                                """
                                                """
                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i]
                                                
                                                
                                            
                                                sc = [item[0] for item in list(seq_score.values())]
                                                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                                max_value_seq = max(max_value_seqs, key=len)
                                                max_strt_idx = seq_score[max_value_seq][1]
                                                """
                                                break
                                            elif score >= S:
                                                temp_seq = seq[i] + temp_seq
                                                score_list.append(score)
                                                seq_score[temp_seq] = [score, i]
                                                break
                                    if backtrack or (backtrack and i == 0) or (not backtrack and (i == len(seq) - 1)):
                                        sc = [item[0] for item in list(seq_score.values())]
                                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                        max_value_seq = min(max_value_seqs, key=len)
                                        max_strt_idx = seq_score[max_value_seq][1]
                                        temp_seq = max_value_seq
                                        score = max(sc) 

                                    if seq_score[temp_seq][0] < S:
                                        continue
                                    else:
                                        
                                        sc = [item[0] for item in list(seq_score.values())]
                                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                                        max_value_seq = min(max_value_seqs, key=len)
                                        temp_seq = max_value_seq
                                        #tar_idx = min(strt_indx, strt_b_indx)
                                        #res = [min(strt_indx, strt_b_indx) + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0])]
                                        if seq_score[max_value_seq][2]:
                                            tar_idx = seq_score[max_value_seq][1]
                                        else:
                                            tar_idx = seq_score[max_value_seq][1] - len(max_value_seq) + 1
                                        if (tar_idx + diff_indx) >= 0:
                                            res = (tar_idx + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0]))
                                            target = seq[:tar_idx] + max_value_seq + seq[tar_idx+len(max_value_seq):]
                                            #target[tar_idx:(tar_idx+len(max_value_seq))] = max_value_seq
                                            if target not in result:
                                                result[target] = [res]
                                                mod_res = [res, temp_seq]
                                                res_seqs[target] = [mod_res]
                                            else:
                                                present = False
                                                res_seq = result[target]
                                                for res_list in res_seq:
                                                    indx = res_list
                                                    if (res[1] + res[2] == indx[1] + indx[2]) and (res[0] + res[2] == indx[0] + indx[2]):
                                                        present = True
                                                        if res[3] >= indx[3]:
                                                            if res[0] < indx[0] and indx in result[target]:
                                                                result[target].remove(indx)
                                                            #result[target].remove(indx)
                                                            if res not in result[target]:
                                                                result[target].append(res)
                                                if not present:
                                                    old_list = result[target]
                                                    if res not in old_list:
                                                        result[target].append(res)
                                                        mod_res = [res, temp_seq]
                                                        res_seqs[target].append(mod_res)                       
        """
        d = dict()
        d['SEQWENCE'] = [(1, 2, 4, 13)]
        """
        
        return result

    def search_two_hit(self, blast_db, *, query=None, pssm=None, T=11, X=5, S=30, A=40):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.
        :param A: max distance A between two hits for the two-hit method.

        :return: dictionary of target sequences and list of HSP tuples.
        """
        t = ttime()
        words = self.get_words(sequence=query, pssm=pssm, T=T)
        """
        for word1 in words:
            for word2 in words:
                match_seqs = self.get_seq_two_hit(word1, word2)
                for target_seq in match_seqs:
        """
        
        
        for word in words:
            seqs = blast_db.get_seq_two_hit(word, self.word_pos)
        
        p = ttime()
        print(t-p)
        result = {}
        #blast_db.match_seqs = json.load(open("/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise5-ge73tag/match_seqs.json"))
        b = len(blast_db.match_seqs)
        for key in blast_db.match_seqs:
            self.hsps = []
            self.checked = []
            #############################
            #if key == "MFSFVDLRLLLLLGATALLTHGQEDIPEVSCIHNGLRVPNGETWKPDVCLICICHNGTAVCDGVLCKEDLDCPNPQKREGECCPFCPEEYVSPDAEVIGVEGPKGDPGPQGPRGPVGPPGQDGIPGQPGLPGPPGPPGPPGPPGLGGNFASQMSYGYDEKSAGVSVPGPMGPSGPRGLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPRGPPGPPGKNGDDGEAGKPGRPGERGPPGPQGARGLPGTAGLPGMKGHRGFSGLDGAKGDTGPAGPKGEPGSPGENGAPGQMGPRGLPGERGRPGPPGSAGARGNDGAVGAAGPPGPTGPTGPPGFPGAAGAKGEAGPQGARGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGPSGPQGPSGAPGPKGNSGEPGAPGNKGDTGAKGEPGPAGVQGPPGPAGEEGKRGARGEPGPSGLPGPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGLTGSPGSPGPDGKTGPPGPAGQDGRPGPAGPPGARGQAGVMGFPGPKGTAGEPGKAGERGVPGPPGAVGPAGKDGEAGAQGAPGPAGPAGERGEQGPAGSPGFQGLPGPAGPPGEAGKPGEQGVPGDLGAPGPSGARGERGFPGERGVQGPPGPAGPRGNNGAPGNDGAKGDTGAPGAPGSQGAPGLQGMPGERGAAGLPGPKGDRGDAGPKGADGSPGKDGVRGLTGPIGPPGPAGAPGDKGETGPSGPAGPTGARGAPGDRGEPGPPGPAGFAGPPGADGQPGAKGEPGDTGVKGDAGPPGPAGPAGPPGPIGNVGAPGPKGSRGAAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGPVGKEGGKGPRGETGPAGRPGEVGPPGPPGPAGEKGSPGADGPAGSPGTPGPQGIAGQRGVVGLPGQRGERGFPGLPGPSGEPGKQGPSGASGERGPPGPMGPPGLAGPPGESGREGSPGAEGSPGRDGAPGAKGDRGETGPAGPPGAPGAPGAPGPVGPAGKNGDRGETGPAGPAGPIGPAGARGPAGPQGPRGDKGETGEQGDRGIKGHRGFSGLQGPPGSPGSPGEQGPSGASGPAGPRGPPGSAGSPGKDGLNGLPGPIGPPGPRGRTGDSGPAGPPGPPGPPGPPGPPSGGYDFSFLPQPPQEKSQDGGRYYRADDANVVRDRDLEVDTTLKSLSQQIENIRSPEGSRKNPARTCRDLKMCHSDWKSGEYWIDPNQGCNLDAIKVYCNMETGQTCVFPTQPSVPQKNWYISPNPKEKKHVWFGESMTDGFQFEYGSEGSDPADVAIQLTFLRLMSTEASQNITYHCKNSVAYMDQQTGNLKKSLLLQGSNEIELRGEGNSRFTYSTLVDGCTSHTGTWGKTVIEYKTTKTSRLPIIDVAPLDIGAPDQEFGMDIGPACFV":
            if True:
            #############################
                seq_list = blast_db.match_seqs[key]
                seq_list_len = len(seq_list)
                for i in range(seq_list_len):
                    indx_list_l = seq_list[i]
                    hl = indx_list_l[2]
                    for j in range(seq_list_len):
                        
                        indx_list_r = seq_list[j]
                        for l in indx_list_l[0]:
                            for r in indx_list_r[0]:
                                for l_seq in indx_list_l[1]:
                                    for r_seq in indx_list_r[1]:
                                        proceed = True
                                        for check in self.checked:
                                            if l == check[0][0] and r == check[1][0]:
                                                proceed = False
                                        ###########

                                        #if l == 100 and r == 104:
                                        if proceed:
                                        ##########
                                        #if (indx_list_l != indx_list_r) and ((indx_list_r[0][r] - indx_list_l[0][l]) in range(3, A+1)) and indx_list_r[0][r] - indx_list_l[0][l] == indx_list_r[1][0] - indx_list_l[1][0]:
                                            if (l != r) and ((r - l) in range(3, A+1)) and (r - l) == (r_seq - l_seq):
                                                #hr = indx_list_r[j]
                                                target, res = self.check_valid_hsp(key, query, [l, l_seq, indx_list_l[2]], [r, r_seq, indx_list_r[2]], X, S, result)
                                                if target != 0 and res != 0:
                                                    if target not in result:
                                                        result[target] = [res]
                                                    else:
                                                        if res not in result[target]:
                                                            found = False
                                                            for res_list in result[target]:
                                                                if res[0] == res_list[0] and res[1] == res_list[1]:
                                                                    found = True
                                                                    if res[2] < res_list[2]:
                                                                        result[target].append(res)
                                                                        result[target].remove(res_list)
                                                                    break
                                                               # elif res[0] in range(res_list[0], res_list[0] + res_list[2] + 1) and (res[0] + res[2]) in range(res_list[0], res_list[0] + res_list[2] + 1) \
                                                                #    and res[1] in range(res_list[1], res_list[1] + res_list[2] + 1) and (res[1] + res[2]) in range(res_list[1], res_list[1] + res_list[2] + 1):
                                                                elif (res[0] + res[2] == res_list[2] + res_list[0])   and (res[1] + res[2]) == (res_list[1] + res_list[2]) :
                                                                
                                                                    found = True
                                                                    break
                                                            if not found:
                                                                result[target].append(res)


        d = dict()
        d['SEQWENCE'] = [(1, 2, 4, 13)]
        p = ttime()
        print(t-p)
        return result
    
    def check_valid_hsp(self, key, query, indx_list_l, indx_list_r, X, S, result):
        strt_indx = indx_list_r[1] + 2
        #strt_indx = indx_list_r[1] 
        end_indx = indx_list_l[1] - 1
        diff_indx = indx_list_r[0] - indx_list_r[1]
        score = 0
        score_list = [0]
        seq_score = {}
        temp_seq = ""
        seq = key
        temp_seq = indx_list_r[2]
        res = 0
        
        ##########
        score = self.sub_matrix[AA_TO_INT[query[strt_indx + diff_indx]]][AA_TO_INT[seq[strt_indx]]] + \
                                        self.sub_matrix[AA_TO_INT[query[strt_indx - 1 + diff_indx]]][AA_TO_INT[seq[strt_indx - 1]]] + \
                                        self.sub_matrix[AA_TO_INT[query[strt_indx - 2 + diff_indx]]][AA_TO_INT[seq[strt_indx - 2]]]
        #print(score)
        score_list.append(score)
        seq_score[temp_seq] = [score, strt_indx - 2]

        #########
        ###########
        #seq = "MGSTLGCHRSIPRDPSDLSHSRKFSAACNFSNILVNQERLNINTATEEELMTLPGVTRAVARSIVEYREYIGGFKKVEDLALVSGVGATKLEQVKFEICVSSKGNSAQHSPSSLRRDLLAEQQPHHLATTVPLTPRVNINTATLAQLMSVRGLSEKMAVSIVDYRREHGPFRSVEDLVRMDGINAAFLDRIRHQVFAERSRPPSTHTNGGLTFTAKPHPSPTSLSLQSEDLDLPPGGPTQIISMRPSVEAFGGMRDGRPVFRLATWNLQGCSVEKANNPGVREVVCMTLLENSIKLLAVQELLDKEALEKFCTELNQPILPNIRKWKGPRGCWRSIVAEKPSSQLQKGPCYSGFLWDTAANVELRDIPGQESSPSNGHAKTVGPSPFLARFKVGSNDLTLVNLQLTALALPGVENSSKNHSDGHRLLNFALTLQETLKGEKDVVILGDFGQGPDSSDYDILRREKFHHLIPAHTFTNISTRNPQGSKSVDNIWISKSLKKVFTGHWAVVREGLTNPWIPDNWSWGGVASEHCPVLAELYMEKDWSKKEVPRNGNGVTLEPSEANVKHER"

        #for i in range(strt_indx - 1, -1 , -1):
        for i in range(strt_indx - 3, -1 , -1):
            score += self.sub_matrix[AA_TO_INT[query[i + diff_indx]]][AA_TO_INT[seq[i]]]                        
            #if ((max(score_list)) - score) < X and (i + diff_indx) != 0 and score < S:
            ##########
            if ((max(score_list)) - score) < X and (i + diff_indx) != 0:
            ##########
                temp_seq = seq[i] + temp_seq
                score_list.append(score)
                seq_score[temp_seq] = [score, i, True]
                strt_b_indx = i
            elif (i + diff_indx) == 0 or (max(score_list) - score) >= X:

                temp_seq = seq[i] + temp_seq
                score_list.append(score)
                seq_score[temp_seq] = [score, i, True]
                break
        #if i in range(end_indx+1, end_indx+4):
        #if indx_list_l[1] in range(i, i+len(temp_seq)):
        if indx_list_l[1]+2 >= i:
            sc = [item[0] for item in list(seq_score.values())]
            max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
            ##########
            #max_value_seq = max(max_value_seqs, key=len)
            max_value_seq = min(max_value_seqs, key=len)
            ##########
            max_strt_idx = seq_score[max_value_seq][1]
            temp_seq = max_value_seq
            score = max(sc) 
            strt_indx = max_strt_idx
            for i in range(strt_indx + len(temp_seq), len(seq)):
                if (i + diff_indx) < len(query):
                    score += self.sub_matrix[AA_TO_INT[query[i + diff_indx]]][AA_TO_INT[seq[i]]]                        
                    #if (score - (max(score_list))) < X and (i + diff_indx) != (len(query) - 1):
                    #if ((max(score_list)) - score) < X and (i + diff_indx) != (len(query) - 1) and score < S:
                    if ((max(score_list)) - score) < X and (i + diff_indx) != (len(query) - 1):
                        temp_seq += seq[i]
                        score_list.append(score)
                        seq_score[temp_seq] = [score, i, False]
                    elif (i + diff_indx) == (len(query) - 1) or (max(score_list) - score) >= X:
                        temp_seq += seq[i]
                        score_list.append(score)
                        seq_score[temp_seq] = [score, i, False]

                        sc = [item[0] for item in list(seq_score.values())]
                        max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                        
                        ##################
                        #max_value_seq = max(max_value_seqs, key=len)
                        max_value_seq = min(max_value_seqs, key=len)
                        ##################
                        max_strt_idx = seq_score[max_value_seq][1]
                        temp_seq = max_value_seq
                        score = max(sc) 
                        break

        ###############
            if (i == len(seq) - 1):
                sc = [item[0] for item in list(seq_score.values())]
                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                #############
                #max_value_seq = max(max_value_seqs, key=len)
                max_value_seq = min(max_value_seqs, key=len)
                ############
                max_strt_idx = seq_score[max_value_seq][1]
                temp_seq = max_value_seq
                score = max(sc) 
            if seq_score[temp_seq][0] < S:
                tar_idx = seq_score[max_value_seq][1]
                self.hsps.append([(tar_idx + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0])), indx_list_r]) 
                self.checked.append([indx_list_l, indx_list_r])
                return 0, 0
            else:
                
                sc = [item[0] for item in list(seq_score.values())]
                max_value_seqs = [key for m in [max(sc)] for key,val in seq_score.items() if val[0] == m]
                
                ###########
                max_value_seq = min(max_value_seqs, key=len)
                #max_value_seq = max(max_value_seqs, key=len)
                ################
                temp_seq = max_value_seq
                #tar_idx = min(strt_indx, strt_b_indx)
                #res = [min(strt_indx, strt_b_indx) + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0])]
                if seq_score[max_value_seq][2]:
                    tar_idx = seq_score[max_value_seq][1]
                else:
                    tar_idx = seq_score[max_value_seq][1] - len(max_value_seq) + 1
                if (tar_idx + diff_indx) >= 0:
                    ##############
                    """
                    if (tar_idx + len(temp_seq) - 1) < indx_list_r[1] + 2
                    mod_seq = seq[tar_idx: ]
                    res = (tar_idx + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0]))
                    """
                    ##############
                    res = (tar_idx + diff_indx, tar_idx,  len(temp_seq), int(seq_score[temp_seq][0]))
                """
                for hsp in self.hsps:
                    #if res[0] in range(hsp[0], hsp[0] + hsp[2]) and (res[0] + res[2]) in range(hsp[0], hsp[0] + hsp[2]) \
                     #   and res[1] in range(hsp[1], hsp[1] + hsp[2]) and (res[1] + res[2]) in range(hsp[1], hsp[1] + hsp[2]):
                    if res[0] == hsp[0][0] and res[1] == hsp[0][1]:
                        hsp_indx_r = hsp[1]
                        if hsp_indx_r[0] < indx_list_r[0]:
                            return 0, 0
                """
                target = seq
        else:
            return 0, 0

        return target, res


"""
blast_db = BlastDB()
from pathlib import Path
import json

test_json = 'blast_test.json'
relative_path = "/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise5-ge73tag/tests/"

with Path(relative_path, test_json).open('r') as json_file:
    json_data = json.load(json_file)
db_sequences =   json_data['db_sequences']
for s in db_sequences:
    blast_db.add_sequence(s)
#blast_db.get_db_stats()
query_seq = json_data['query_seq']
blast = Blast(json_data['sub_matrix'])
blast_ex = json_data["blast_hsp_one_hit"]
#pssm = json_data["query_pssm"]
#words = blast_db.get_words(sequence=query, T=11)
results = blast.search_two_hit(blast_db,
                                       query=query_seq,
                                       pssm = None,
                                       T=11,
                                       X=5,
                                       S=30, A=40)
"""