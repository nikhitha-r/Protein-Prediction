from itertools import cycle

def get_orfs(genome, min_num_aa):
    """
    Find and return all ORFs within the provided genome sequence.

    :param genome: String with the genome sequence.
    :param min_num_aa: Minimum number of amino acids in the translated protein
                       sequence.

    :return: List of tuples representing the ORFs (details in worksheet).
    """
    allowed = "ATGC"
    if not isinstance(genome, str) or len(genome) == 0 or not all(c in allowed for c in genome):
        raise TypeError
    start_codon = "ATG"
    stop_codon = ['TAA', 'TAG', 'TGA']
    ref_dict = {"T" : "A", "A" : "T", "G" : "C", "C" : "G"}
    amino_dict = {
            'L' : ['CTC', 'CTT', 'CTA', 'CTG', 'TTA', 'TTG'],
            'S' : ['TCA', 'TCT', 'TCC', 'TCG', 'AGC', 'AGT'],
            'R' : ['CGA', 'CGC', 'CGT', 'CGG', 'AGA', 'AGG'],
            'V' : ['GTA', 'GTG', 'GTC', 'GTT'],
            'P' : ['CCC', 'CCA', 'CCG', 'CCT'],
            'T' : ['ACC', 'ACG', 'ACT', 'ACA'],
            'A' : ['GCA', 'GCC', 'GCG', 'GCT'],
            'G' : ['GGA', 'GGC', 'GGT', 'GGG'],
            'I' : ['ATA', 'ATC', 'ATT'],
            'F' : ['TTT', 'TTC'],
            'Y' : ['TAT', 'TAC'],
            'H' : ['CAC', 'CAT'],
            'Q' : ['CAG', 'CAA'],
            'N' : ['AAC', 'AAT'],
            'K' : ['AAA', 'AAG'],
            'D' : ['GAC', 'GAT'],
            'E' : ['GAA', 'GAG'],
            'C' : ['TGC', 'TGT'],
            'M' : ['ATG'],
            'W' : ['TGG']

        }
    comp_genome = ""
    for stra in genome:
        comp_genome += ref_dict[stra]
    main_orfs = find_orfs(genome, start_codon, stop_codon, min_num_aa, amino_dict, False)
    comp_orfs = find_orfs(comp_genome[::-1], start_codon, stop_codon, min_num_aa, amino_dict, True)
    circular_orfs = find_cir_orfs(genome, main_orfs, start_codon, stop_codon, min_num_aa, amino_dict, False)
    
    circular_orfs_comp = find_cir_orfs(comp_genome[::-1], comp_orfs, start_codon, stop_codon, min_num_aa, amino_dict, True)
    
    for main_orf in main_orfs:
        for cir_orf in circular_orfs:
            if main_orf[0] <= cir_orf[1] and main_orf[1] <= cir_orf[1] or len(main_orf) == 5:
                main_orfs.remove(main_orf)
    for comp_orf in comp_orfs:
        for cir_orf in circular_orfs_comp:
            if  comp_orf[1] == cir_orf[1] or len(comp_orf) == 5:
                comp_orfs.remove(comp_orf)

    final_orf = main_orfs + comp_orfs + circular_orfs + circular_orfs_comp
    #print(len(comp_orfs))
    
    
    
    return final_orf

def find_cir_orfs(genome, main_orfs, start_codon, stop_codon, min_num_aa, amino_dict, isReverse):
    all_circ_ofs = []
    for orf in main_orfs:
        if len(orf) == 5:
            start_idx = orf[0]
            mod_list = genome[start_idx:] + genome[:start_idx]
            thres = len(genome) - start_idx
            cir_orfs = find_orfs(mod_list, start_codon, stop_codon, min_num_aa, amino_dict, isReverse)
            #cir_orfs = find_orfs(mod_list, "GTA", ['AAT', 'GAT', 'AGT'], min_num_aa, amino_dict, False)
            for cir_orf in cir_orfs:
                if isReverse:
                    indx = (cir_orf[0]+thres - 1)%(len(genome) -1 ) 
                    if indx < thres:
                        cir_orf_mod = (indx, cir_orf[1] +  thres, cir_orf[2], cir_orf[3])
                        all_circ_ofs.append(cir_orf_mod)
                elif not isReverse and cir_orf[0] <= thres:
                    cir_orf_mod = (start_idx+cir_orf[0], (cir_orf[1] +  start_idx+cir_orf[0] - 1)%(len(genome) -1 ), cir_orf[2], cir_orf[3])
                    all_circ_ofs.append(cir_orf_mod)
    return all_circ_ofs




def find_orfs(genome, start_codon, stop_codon, min_num_aa, amino_dict, isReverse):
    orfs = []
    gen_length = len(genome)
    for i in range(3):
        orf = ""
        aa_len = 0
        start_idx = None
        sub_seq_track = {}
        """
        if isReverse:
            start = gen_length - 1
            end = i + 1
            step = -3

        else:
            """ 
        start = i
        end = gen_length
        step = 3
        for j in range(start, end, step):
           # if isReverse:
           #     sub_seq = genome[j + step:j]
           # else:
            sub_seq = genome[j:j+ step]
            if (sub_seq == start_codon) or (len(orf) != 0) and (sub_seq not in stop_codon):
                if j < end-1 and j > end-4 and isReverse:
                    isCircular = True
                    orf_list = (start_idx, end_idx, orf, False, isCircular)
                    orfs.append(orf_list)
                if j < end-1 and j > end-4 and not isReverse:
                    isCircular = True
                    orf_list = (start_idx, end_idx, orf, False, isCircular)
                    orfs.append(orf_list)
                if (sub_seq == start_codon) and (start_idx is None):
                    start_idx = j
                aa_len += 1
                orf += find_amino_acid(sub_seq, amino_dict)
            elif sub_seq in stop_codon:
                
                end_idx = j + (step - 1)
                
                #orf += sub_seq
                if aa_len >= min_num_aa:
                    
                    if isReverse:
                        orf_list = (gen_length - start_idx - 1, gen_length - end_idx - 1, orf, True)
                    else:
                        orf_list = (start_idx, end_idx, orf, False)
                    """
                    if start_idx not in sub_seq_track:
                        sub_seq_track[start_idx] = [sub_seq]               
                        orfs.append(orf_list)
                        ###########
                        for li in orfs:
                            if li[0] == start_idx and (li[2][-3:] == sub_seq) and (abs(li[1] - li[0]) < abs(start_idx - end_idx)) and len(li) == 4:
                                orfs.remove(li)
                                orfs.append(orf_list)
                            elif li[0] > start_idx and li[1] < end_idx and len(li) == 4:
                                orfs.remove(li)
                                orfs.append(orf_list)
                            else:
                                orfs.append(orf_list)
                        ##############
                    else:
                        stop_seq_list = sub_seq_track[start_idx]
                        if sub_seq in stop_seq_list:

                            for li in orfs:
                                if li[0] == start_idx and (li[2][-3:] == sub_seq) and (abs(li[1] - li[0]) < abs(start_idx - end_idx)) and len(li) == 4:
                                    orfs.remove(li)
                                    orfs.append(orf_list)
                                elif li[0] > start_idx and li[1] < end_idx and len(li) == 4:
                                    orfs.remove(li)
                                    orfs.append(orf_list)
                                else:
                                    orfs.append(orf_list)
                        else:
                            sub_seq_track[start_idx].append(sub_seq)
                            orfs.append(orf_list)
                        
                    """
                    if start_idx not in sub_seq_track:
                        sub_seq_track[start_idx] = [sub_seq]  
                    else:
                        sub_seq_track[start_idx].append(sub_seq)
                    isAdded = False
                    for st_idx in sub_seq_track:
                        if st_idx >= start_idx and (sub_seq in sub_seq_track[st_idx]) and not isReverse:
                            for li in orfs:
                                if li[0] == st_idx and li[1] < end_idx and (abs(li[1] - li[0]) < abs(start_idx - end_idx)) and len(li) == 4:
                                    orfs.remove(li)
                                    orfs.append(orf_list)
                                    isAdded = True

                        elif st_idx <= start_idx and (sub_seq in sub_seq_track[st_idx]) and isReverse:
                            for li in orfs:
                                if li[0] == st_idx and li[1] < end_idx and (abs(li[1] - li[0]) < abs(start_idx - end_idx)) and len(li) == 4:
                                    orfs.remove(li)
                                    orfs.append(orf_list)
                                    isAdded = True   
                    if not isAdded:
                        orfs.append(orf_list)
                #break
                orf = ""
                aa_len = 0
                start_idx = None
    return orfs

def find_amino_acid(seq, amino_dict):
    
    for key in amino_dict:
        seq_list = amino_dict[key]
        if seq in seq_list:
            return str(key)
    return ""

get_orfs("AGATGCATCTAAATACCAGAATGGATTATCTGTTTGTTTCGTAATCCACTGATCTAGCGCCGAACGAACGCCTATCTCAACTTTTTCAATTATATCAAAGATATAACTTCGCAACTCACTATCAAACATATAAAGTTCATAACCATTATGAAACGTAGTATGTTCTGAAAATTTCTTAGTCCCTTTATCTTTAAAGGGTGATAGATAAGCTTTAAATCGATAGTAACTACAACGATTTAAGACTTTTTCAGCAAAATCTTCATCATCAATGATTAGCCCTTGGTCTATGAGTTTTTTACATAACTGAGAACTAGATAGGTAAGGTTTTTTATAAGGAACGAGGGGAGCCGTAGTCATATCCATTTTCTCTAGGCAAAAAAAACTTCCACTCTTACGAACTTGCGTTCAGGAGGAAGTCATTGTTAAGTTAATTATACATGTTACATATACAGAGTCCAGTTATTTTCAATCACGTGACTCCGTAACGCAGAACCTAAGATCCCCGCGGGTGAAGTAGGTTGTCTTGTGTTGAATCTAGGTACGCACAAAATTGGTGTTCATTATGCCCCTTATGTTAAATCCGTAAGTTAGCCTATAAAAACCATGCCTTCGCAGCCTCAGAAGCCCATAGAGGCGTTTTTAAGATATAAACCACACCAACATACCCCAAAACCAGTTTGAAAGTTTATGACGCTTTCTAGGCGATTTACAGGGCAATAAAATCACAATCCGTTTCAATCAATAATCACCCTGTAACGTCTATTTTCAACAATTAGCCGGTGCGCTCGAAAAAGCCCGAAGCGTAGCTGGCTTTTTTGGGCGTACTGGAGACAAACGAAGTGCGTCAGGGTTTCAAAGGTTTCACCTTTGGCACGGGATGTTTTTATCCAATGCGTAGCAGGGGGTGAAAACGTCCCCGTCGTGCCTATCTTGAAAATAGAAATCAACGGGTAATCACAACGTAACTAACATTAAAAAGCTAACCCATTGAATCAATGAGTTTATAGGTAGGTTACTAAGATAAATCATGCCAATTTCAGCCCAACTGAAATAACCATTCTCGTTCTTCTTTTTTATTGCTAGTCTCTAAAGTATCAATCAGTAATCAATAAAGTATTCAGTGAAAAAGAAAACACCGTTAGAAAAAGAGCTTCTTTCGTACACGAAAAAGCATCAAACGAAATCTAAAACCGCCGTCATTACGGAGTTATATGACGTGATAACCATCTTGATGAGTGAAGGTATCTCACTAGAGCACATCGTTGCGTTTTTAAAGACAAAAAATATCACGGTATCGAGAGACTATTTGAAGACGGTTCTTTATCGAATTAGGCAAAAGAGAAAACTAGAAAATGAATAAACCAAAATCCAAGACTCATTTGGTGTCCGTTCGAGTGGATGGTCGAGTGAATGCTGAAATTGAAGAGCTGTTAAACGACACGCATAAAGATAAAGCAGAACTGTTAAGGGCGTTGCTAGTAAGAGGTTTGTTCGCCTATAAAAGGGATTTGTATCAAACTCATGGCGATAGAGTGCATCACTCAATTGATAAATCGTATTACGCCTACGAAAGAGAAATGCTTGATGATGACTTACATTACATCGCCACTCGTTTTGACCAACTAGAAGCCATGATAAGAGAGCATAAAGAGCCGGTATTAACACCAACGAAAAAAGAAGGAAAAAAGACGCTTTTGGATAGTTTCTTTAAAAAATAACGTTACCCATAAATTAACGTATTAATACAAGGAGTTAGCATTTTATTATTGATTACTCAGTGATATCTATTTTTTAGATAGGCACGACGGGCACGTTTTCACCCTCTACATACGTATTTGGAAAAACATGCCGCGCCAAAGGATGTCCTTTGAAACCTAAGAGCAAAAAAGAGAATGAAGAAAAGACGATTTACTCAAGAAGACCGCACTCTTATTGAACTACAAATAGAAAAAGGAGTGGCACGAAAAGAAGCGGCCAAAGCAATTTATGAAAAAGTATACGGACTCAAGCCAAGAGACTTTGCTCTTATAGATTCACTCTGCAGCAGCAACGGTGAAGATTATAGAGACCAAGCCTTTGAAACGATTTTCAAAAGGAACAGAGAAAAATTAAAACGACTCACCATAGCAAGGGCAGCGAGAGCGGAATTGAAAAAGTTACCGAAAGAAGAAAAAAGCACATTCAGTTGGTCAGCATCGATCCATAAAAGAAGATAATCGACCCAACTAACGCGCTTCGCTTGTTCGCACATTAGAAAAGAAAGTGGGCTTAACTCGCTTATGCTCGTTGCGGCAACTTTCATTTCTAATCGTGCAGCTCTTAATATTTTCTTTTTAATTCTTTTTGAATTGAAGAGCAAAAGAAGTCGGCATGTATGCTGATAGGTGAGTTTTGAGCTAGATTTCGGTAATACTCACCCTAGAGGATAACCGCTTTCGCTAAGTTTATCACTCTAGTGTGAACCCTGATTGCTTGCCGACTGGAGCCGTTCATCGCCATCATTCGCTTGAATTCTAATACACTATTCAAGTGCAGAGATTACGGGCTGTATTGCCCCCATCCTCACCGTTCTCTGCCGACCTCTGCATCAGTACCTTGCAAGGGTTTCCCCAATCTACAAGCCTCTTTGCCTAGCTGCTTACGCTAGTTCTTGCTTCAGTCCACTCATTTGGTGTGCGTTATCGGATGTTACGCTAGTCTCTCGTCCAGTATCTGGGCGAAATATGGTATGAAAACCAACCGCATGATGCTGTCCCACAACCCTAGTGACTCTAAACCCAACCACTGATTGCTTGGAAATACTTGCAGTTCTTATGTGATATGATCTATCATGAAAGGACTTACAACGCCGTCTAAACGTAGTTAAGTATTTCGAAAAAACCAAGCTTGCAGGCTTGGTTTTTTCACATCTGGGCATAAGGTTACGGCTAAATGCTCAAGATCTCAATGATCAATTTTATGATCCTCTGGTTCATAGACAAAAAAGCCACTCTTTTGAGTGGCTTTTTTGTTTCTATTCTTTGTGTGATTTTACATACATTAAATGAATACCTGTCGGGGTTATTTGGGTTAGGTTATAACTGTATTCATAACCGTTTTCATCAAGCAATTTTATACAATATAGTGCAAATACAACCGCAAATACTACACCTAATATCATGCTATAAATACACGTTCTATTGATAGTACTAGTCCATTTTTCACCAATGACCAAACGAGCTTGTTTCCCCGAATAAATAGCAATAGAAATTAAAGCAACGGAAATAAGTAAAGGGATAAAACCACCAATGGGGACGATTAATGTTATCAGTGATAATCTAACTTCATCTGTCTGATTTAATAAATCTTGGTAATCCCCCAAGACTCCATTTGTAGGTATCCACCAAAGTAAAATCAATAAAGGTATAAATACTATTCCTAAAAACCTCATTTTAACTGATAACACTTTCTAATCCCTTAACTACCTCATCATTTACCTTGTATTCACTCCATATCCACTCAATAGACACATAGACAAAAAACACAGCAATCAAGCCCAGAGTAACGTAACCAGTCATCGCAATAATTCCAAACAATATAGCCTCTCCTGCACCAAACTGTAGCAATGCCTTAAACATATCTGCCCCAACATTACCGAACCAATCAACTAAGTGATAATCGTCTTTAAAAACAAGATCTGTCGTTGCTATCGCTGCCGATACCACAAACGTCAATACACCTGCGCCTTTAAATCCATTTGCTCGTGCTTTAGGACTCAAACCGACTTGTTGAACTTTAGGGCTGTTAATTGGGTACTTCTTACCGTTCATGTTAATTCGAGTGCCATTAACAAGTGCGTGTTAAATTTCTTTGCCGTTCTTCTTACCAGAGAAAGCAATAAACTGTTTCCCCTTACTCTCTACAATATCCGCTTTAATATTTATATCGCCGAATTGCTTAACTAAACCTTTAGCATCAAATATCGGTGTTGTATGTTGAAACCAACCGTTAGCCGTATTTAAACTTGTTGCTATTAACGCTGTATCTTTGCCTAAATCCACCCAAAAATTATTTAATAATTTATGTGAGTCTTCCAAATCTAAAATAACAACATCTTGGTGATTTTGGTTTAATTCATCCTCTAAAGAGCACTGTAAACCCCAAGAATCATCTAATGCCTGAATCGCTTTAGAACCTTTAGTGACGAGTGGGAAATCATAAGCTTGTTGCGTAACTAATTGAGCTTGAGTTTCGTTTGCTTGAATATGACTGTTTGAACCGCCATTTAACGATTTGGTAAACAAAGAAAGAGGAGATACACACTTTTCACCAAAAGGCTTTCGATATGCAGTTAATGACTCTAAACCATCCTCTGTTGGTACAATACGATTCCGAATATTATCATTGAGAATTTTAGAAGCCTTACACCCTCGATTTTCATCTGCTAAATCCCTTTGCTCTCTCACTTCATACTCTGCATTCATCGTTTCAGTAAACACTTGTTGAGCTTCTTGGAGCGTTGCGTATGTCACACCTCTAAATACTTTTGAATAAACTAAAGGATCATGCCCCAAATTCTGTTGTGTTGCGGTCTGATTAAAAAACATATCAACGGCTTGGTCATACGTTACATGACTCACTAACCAATCAGCAGAGTGACGCGGTGCAAATGGAAAACTTCGACCTTTAGGTGGCGGCATATACAAACTCCAATAAAAAATAATTAAATATTACCAGATAGCTAAAAATCAAAACGCCTGAATAAGCCTGTTAAAAATACTATTTCGTTAGGACAAACGACATTTTTATCACTTAAAAAACAGTTTTCATAAAAAAAGGCCACCTAAATAGGCGACCTTTTACGATTCGTAGAACGGGTTTTATGTTAAGTATAGCGATTATAAATCAAAAAAAAGCGGCTCTTCTTTCCAACCATTAGGGAAACCCATTGATAAGCTAAACCGCTGTGTTGTTGTATAATTATCAAATAAATCAGAGATTTTTGGTCCTATTTTTTCATCGAAACCAAGCCCTGAATATATACATTGTAATGCAGCAAGAGCTGTATAAAGTCGGTTTAACTGATCCTCTTCACGCTTATCTGGATTTGGTTTGGTCCTGACAAGAGGAATTGCATCAGACAAGCTCTTTTTAATCGCTGTTGGAGCAGGGAAGTTTCGATTAAACAGTCGGTTATGATGCCCACAGTAGTTTCTCATCTGGTGCAGCACATTCATCCAACTAATTAATGTCTGAAATTTCTGAATATTGAACTTCTTAGAAAATCTATCCATCTTAAGTGATTGAATTTTATCATCACTGATGTTTTGAATTAACTTCACTACGTTGCCGAAAGTCATTAACTCAATGGCTACCCATGCAGGCGGGAGATCTCTGTAAAACGGACAATACTCATTGTAATACTTTTCTTGATAGTGCTTTGCGTACTCTTCTTTAGAGTCAACAAACATATCCCTAACCCTGCTAACCGTTTTTACTTGTTGGCCATTGCTCTTAAATAA", 34)
