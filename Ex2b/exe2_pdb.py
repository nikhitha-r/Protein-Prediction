# -*- coding: utf-8 -*-

from Bio.PDB.MMCIFParser import MMCIFParser # Tip: This module might be useful for parsing... 
import numpy as np


############# Exercise 2: Protein Data Bank #############
# General remark: In our exercise every structure will have EXACTLY ONE model.
# This is true for nearly all X-Ray structures. NMR structures have several models.
class PDB_Parser:
    
    CIF_PARSER     = MMCIFParser() # parser object for reading in structure in CIF format
    
    def __init__( self, path ):
        '''
            Initialize every PDB_Parser with a path to a structure-file in CIF format.
            An example file is included in the repository (7ahl.cif).
            Tip: Store the parsed structure in an object variable instead of parsing it
            again & again ...
        '''


        self.structure = self.CIF_PARSER.get_structure("pdb", path) # Parse the structure once and re-use it in the functions below
        self.one_letter_code = {'ALA': "A", "ARG": "R", "ASN": "N", "ASP":"D", "ASX":"B", "CYS": "C", "GLU" : "E",
                                "GLN" : "Q", "GLX" :"Z", "GLY": "G", "HIS": "H", "ILE": "I", "LEU" : "L",
                                "LYS" : "K", "MET": "M", "PHE": "F", "PRO" : "P", "SER" : "S", "THR": "T",
                                "TRP": "W", "TYR": "Y", "VAL": "V"}
    # 2.8 Chains    
    def get_number_of_chains( self ):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
            Return:
                Number of chains in this structure as integer.
        '''
        n_chains = len(self.structure.child_dict[0].child_dict)
        return n_chains
        
    
    # 2.9 Sequence  
    def get_sequence( self, chain_id ):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the amino acid sequence (single-letter alphabet!) of a given chain (chain_id)
                in a Biopython.PDB structure as a string.
        '''
        chain_dict = self.structure.child_list[0].child_dict[chain_id]
        child_list = chain_dict.child_list
        
        sequence = ""
        for res in child_list:
            if res.resname in self.one_letter_code:
                one_let_code = self.one_letter_code[res.resname]
            
                sequence += one_let_code
        #sequence = self.CIF_PARSER._mmcif_dict['_entity_poly.pdbx_seq_one_letter_code'][0].replace('\n', "")
        return sequence
        
        
    # 2.10 Water molecules
    def get_number_of_water_molecules( self, chain_id ):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the number of water molecules of a given chain (chain_id)
                in a Biopython.PDB structure as an integer.
        '''
        chain_dict = self.structure.child_dict[0].child_dict[chain_id]
        res_dict = chain_dict.child_dict
        n_waters = 0
        for key in res_dict.keys():
            if key[0] == 'W':
                n_waters += 1
        return n_waters
    
    
    # 2.11 C-Alpha distance    
    def get_ca_distance( self, chain_id_1, index_1, chain_id_2, index_2 ):
        ''' 
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id_1 : String (usually in ['A','B', 'C' ...]. The number of chains
                                depends on the specific protein and the resulting structure)
                index_1    : index of a residue in a given chain in a Biopython.PDB structure
                chain_id_2 : String (usually in ['A','B', 'C' ...]. The number of chains
                            depends on the specific protein and the resulting structure)
                index_2    : index of a residue in a given chain in a Biopython.PDB structure
        
                chain_id_1 and index_1 describe precisely one residue in a PDB structure,
                chain_id_2 and index_2 describe the second residue.
        
            Return: 
                Return the C-alpha (!) distance between the two residues, described by 
                chain_id_1/index_1 and chain_id_2/index_2. Round the returned value via int().
            
            The reason for using two different chains as an input is that also the distance
            between residues of different chains can be interesting.
            Different chains in a PDB structure can either occur between two different proteins 
            (Heterodimers) or between different copies of the same protein (Homodimers).
        '''
        res1_coord = self.structure.child_dict[0].child_dict[chain_id_1].child_dict[(' ', index_1, ' ')].child_dict["CA"].coord
        res2_coord = self.structure.child_dict[0].child_dict[chain_id_2].child_dict[(' ', index_2, ' ')].child_dict["CA"].coord
        #res1_coord = self.structure.child_dict[0].child_dict[chain_id_1].child_list[index_1].child_dict["CA"].coord
        #res2_coord = self.structure.child_dict[0].child_dict[chain_id_2].child_list[index_2].child_dict["CA"].coord
        diff = np.array(res1_coord) - np.array(res2_coord)

        ca_distance = np.sqrt(np.sum(diff * diff))
        return int( ca_distance )
        
        
    # 2.12 Contact Map  
    def get_contact_map( self, chain_id ):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return a complete contact map (see description in exercise sheet) 
                for a given chain in a Biopython.PDB structure as numpy array. 
                The values in the matrix describe the c-alpha distance between all residues 
                in a chain of a Biopython.PDB structure.
                Only integer values of the distance have to be given (see below).
        '''
        sequence = self.get_sequence(chain_id)
        length = len(sequence)
        if length == 0:
            return None
        contact_map = np.zeros( (length,length), dtype=np.float32 )
        for seq1 in range(length):
            for seq2 in range(length):
                #print(sequence)
                res1_coord = self.structure.child_dict[0].child_dict[chain_id].child_list[seq1].child_dict["CA"].coord
                res2_coord = self.structure.child_dict[0].child_dict[chain_id].child_list[seq2].child_dict["CA"].coord
                diff = np.array(res1_coord) - np.array(res2_coord)

                ca_distance = np.sqrt(np.sum(diff * diff))
                contact_map[seq1][seq2] = ca_distance
                #contact_map[seq1-1][seq2-1] = self.get_sequence(chain_id)
                #contact_map[seq1-1][seq2-1] = self.get_ca_distance(chain_id, seq1, chain_id, seq2 )

        return contact_map.astype( np.int64 ) # return rounded (integer) values
        
        
    # 2.13 B-Factors    
    def get_bfactors( self, chain_id ):
        '''
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the B-Factors for all residues in a chain of a Biopython.PDB structure.
                The B-Factors describe the mobility of an atom or a residue.
                In a Biopython.PDB structure B-Factors are given for each atom in a residue.
                Calculate the mean B-Factor for a residue by averaging over the B-Factor 
                of all atoms in a residue.
                Sometimes B-Factors are not available for a certain residue; 
                (e.g. the residue was not resolved); insert np.nan for those cases.
            
                Finally normalize your B-Factors using Standard scores (zero mean, unit variance).
                You have to use np.nanmean, np.nanvar etc. if you have nan values in your array.
                The returned data structure has to be a numpy array rounded again to integer.
        '''
        length = len(self.get_sequence(chain_id))
        b_factors = np.array( [])
        res_list = self.structure.child_dict[0].child_dict[chain_id].child_list
        for res in res_list:
            if res.resname != "HOH":
                b_fact_res = []
                atom_list = res.child_list
                for atom in atom_list:
                    b_fact_res.append(atom.bfactor)
                b_factors = np.append(b_factors, np.mean(np.array(b_fact_res)))
        b_factors -= np.mean(b_factors)
        b_factors /= np.std(b_factors)
        return b_factors.astype( np.int64 ) # return rounded (integer) values
            
    
def main():
    print('PDB parser class.')
    return None
    

if __name__ == '__main__':
    main()

"""
p = PDB_Parser("/Users/nikhitha/Documents/Protein Prediction/Exercises/pp1cs2020exercise2-ge73tag/tests/7ahl.cif")
p.get_bfactors("A")
"""

