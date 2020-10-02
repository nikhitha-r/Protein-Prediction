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
    #if not all(isinstance(item, int) for item in sequence) or len(sequence) == 0 or (not isinstance(sequence, list)) or sequence == None:
    #    return None, None
    if sequence is None:
        return None, None
    if not isinstance(sequence, list):
        return None, None
    if len(sequence) == 0:
        return None, None
    if not all(isinstance(item, int) for item in sequence):
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
    f_asc_lists = []
    for max_ind in indices:
        f_asc_list = []
        for i in range(max(len_list) - 1):
            if len(f_asc_list) == 0:
                f_asc_list.append(max_ind)
                max_ind = asc_list[max_ind]
                f_asc_list.append(max_ind)
            else:
                f_asc_list.append(asc_list[max_ind])
                max_ind = asc_list[max_ind]
        f_asc_lists.append(f_asc_list)
    final_asc_lists = []
    for f_li in f_asc_lists:
        final_asc_list = []
        for i in range(len(f_li)-1, -1, -1):
            final_asc_list.append(sequence[f_li[i]])
        final_asc_lists.append(final_asc_list)
    final_asc_lists.sort()
    asc_res = list(k for k,_ in itertools.groupby(final_asc_lists))[0]
    """
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
    """  
    if len(asc_res) == 0:
        asc_res.append(sequence[-1])    
    print(asc_list, asc_res, len_list)


    ## DESC
    len_desc_list = []
    desc_list = []
    for i in range(len(sequence)):
        len_desc_list.append(1)
        desc_list.append(0)
    for i in range(1, len(sequence)):
        for j in range(i):
            if sequence[j] > sequence[i]:
                if (1+len_desc_list[j]) > len_desc_list[i]:
                    len_desc_list[i] = 1 + len_desc_list[j]
                    desc_list[i] = j
                elif (1+len_desc_list[j]) == len_desc_list[i]:
                    desc_list[i] = j
    indices = [i for i, x in enumerate(len_desc_list) if x == max(len_desc_list)]
    f_desc_lists = []
    for max_ind in indices:
        f_desc_list = []
        for i in range(max(len_desc_list) - 1):
            if len(f_desc_list) == 0:
                f_desc_list.append(max_ind)
                max_ind = desc_list[max_ind]
                f_desc_list.append(max_ind)
            else:
                f_desc_list.append(desc_list[max_ind])
                max_ind = desc_list[max_ind]
        f_desc_lists.append(f_desc_list)
    final_desc_lists = []
    for f_li in f_desc_lists:
        final_desc_list = []
        for i in range(len(f_li)-1, -1, -1):
            final_desc_list.append(sequence[f_li[i]])
        final_desc_lists.append(final_desc_list)
    final_desc_lists.sort()
    desc_res = list(k for k,_ in itertools.groupby(final_desc_lists))[-1]
    if len(desc_res) == 0:
        desc_res.append(sequence[-1])
    """
    max_len_indx = max(indices)  
    f_desc_list = []
    for i in range(max(len_desc_list) - 1):
        if len(f_desc_list) == 0:
            f_desc_list.append(max_len_indx)
            max_len_indx = desc_list[max_len_indx]
            f_desc_list.append(max_len_indx)
        else:
            f_desc_list.append(desc_list[max_len_indx])
            max_len_indx = desc_list[max_len_indx]
    final_desc_list = []
    for i in range(len(f_desc_list)-1, -1, -1):
        final_desc_list.append(sequence[f_desc_list[i]])
     """       
    print(desc_list, desc_res, len_desc_list)
    return asc_res, desc_res
    #return [], []
"""
print(get_longest_subsequences({1:[-5,
            -4,
            -3,
            -2,
            -1,
            0,
            1,
            2,
            3,
            4]}))
  """       
"""get_longest_subsequences([0, 4, 12, 2 ,10,6, 9, 13, 3, 11, 7, 15])"""