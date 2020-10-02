# This function extracts the longest increasing and decreasing subsequences from
# a sequence of integers.
# The function takes a list of integers and returns two lists of integers 
# containing the longest subsequences. If it gets passed an empty list, None,
# or something other then an integer list, the function returns two times None.
#   
# The subsequence length is not necessarily unique, in this case the elements of
# the longest increasing subsequence should be as small as possible, while the 
# elements of the longest decreasing subsequence should be as large as possible.

import itertools
def get_longest_subsequences(sequence):
    if not all(isinstance(item, int) for item in sequence) or len(sequence) == 0 or (not isinstance(sequence, list)):
        return None, None
    len_list = []
    asc_list = []
    for i in range(len(sequence)):
        len_list.append(1)
        asc_list.append(0)
    for i in range(1, len(sequence)):
        for j in range(i):
            if sequence[j] < sequence[i]:
                if (1+len_list[j]) > len_list[i]:
                    len_list[i] = 1 + len_list[j]
                    asc_list[i] = j
                elif (1+len_list[j]) == len_list[i]:
                    asc_list[i] = j
    indices = [i for i, x in enumerate(len_list) if x == max(len_list)]

    max_len_indx = max(indices)  
    f_asc_list = []
    for i in range(max(len_list) - 1):
        if len(f_asc_list) == 0:
            f_asc_list.append(max_len_indx)
            max_len_indx = asc_list[max_len_indx]
            f_asc_list.append(max_len_indx)
        else:
            f_asc_list.append(asc_list[max_len_indx])
            max_len_indx = asc_list[max_len_indx]
    final_asc_list = []
    for i in range(len(f_asc_list)-1, -1, -1):
        final_asc_list.append(sequence[f_asc_list[i]])
            
    print(asc_list, final_asc_list, len_list)
    return [], []

get_longest_subsequences([-1599,
            4755,
            -9887,
            9436,
            -5354,
            -9141,
            659,
            4247,
            -5844,
            -1521])
            
"""get_longest_subsequences([0, 4, 12, 2 ,10,6, 9, 13, 3, 11, 7, 15])"""