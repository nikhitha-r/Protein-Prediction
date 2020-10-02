# The following function should take a string and a sub-string and return 
# a list of starting positions for the substring

def get_sub_position(search_string, needle):
    search_len = len(needle)
    
    
    indx_list = []
    for j in range(len(search_string)):
        ne_char = search_string[j]
        if ne_char == needle[0]:
            if needle == search_string[j: j+search_len]:
                indx_list.append(j)


    return indx_list


