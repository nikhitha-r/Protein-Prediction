# The following method should take a strand as an argument and
# return a complementary strand

def complementary(strand):
    ref_dict = {"T" : "A", "A" : "T", "G" : "C", "C" : "G"}
    str_list = ""
    for stra in strand:
        str_list += ref_dict[stra]
    return str_list
    
    