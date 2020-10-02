# -*- coding: utf-8 -*-
from Bio import SeqIO # Tip: This module might be useful for parsing... 

############ Exercise 3: SwissProt ##########
class SwissProt_Parser:
    
    PARSER = SeqIO
    
    def __init__( self, path, frmt='uniprot-xml' ):
        '''
            Initialize every SwissProt_Parser with a path to a XML-formatted UniProt file.
            An example file is included in the repository (P09616.xml).
            Tip: Store the parsed XML entry in an object variable instead of parsing it
            again & again ...
        '''
        self.sp_anno = SeqIO.parse(path, frmt) # Parse the XML file once and re-use it in the functions below
        for record in self.sp_anno:
            self.sp_record = record 

    # 2.2 SwissProt Identifiers
    def get_sp_identifier( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Unique SwissProt identifier for the given xml file
        '''
        identifier = self.sp_record.id
        return identifier
        
    # 2.3 SwissProt Sequence length
    def get_sp_sequence_length( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Return sequence length of the UniProt entry as an integer.
        '''
        #print("SADDS12121212121", self.sp_anno)
        seq = self.sp_record.seq
        seq_len = len(seq)
        
        return seq_len
    
    # 2.4 Organism 
    def get_organism( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Return the name of the organsim as stated in the corresponding field
                of the XML data. Return value has to be a string.
        '''
        organism = self.sp_record.annotations['organism']
        return organism
    
    # 2.5 Localizations
    def get_localization( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Return the name of the subcellular localization as stated in the 
                corresponding field.
                Return value has to be a list of strings.
        '''
        localization = self.sp_record.annotations["comment_subcellularlocation_location"]
        return localization
    
    # 2.6 Cross-references to PDB
    def get_pdb_support( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Returns a list of all PDB IDs which support the annotation of the
                given SwissProt XML file. Return the PDB IDs as list.
        '''
        dbrefs = self.sp_record.dbxrefs
        pdb_ids = []
        for ref in dbrefs:
            if 'PDB:' in ref:
                pdb_ids.append(ref[4:])            
        return pdb_ids
        
def main():
    print('SwissProt XML Parser class')
    return None
    
if __name__ == '__main__':
    main()
"""
s = SwissProt_Parser("/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise2-ge73tag/tests/P09616.xml")
s.get_sp_identifier()
s.get_sp_sequence_length()
"""

