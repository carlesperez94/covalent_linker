import os 

from rdkit import Chem
from covalent_linker.Reactions.pdb import PDB, get_specifyc_atom_pdb_name, get_resnum_from_line
from covalent_linker.Reactions.pdb import get_atom_pdb_name_from_line, get_resname_from_line
from covalent_linker.Helpers.filemanager import create_file 

class Reactant:
    """
    It defines reactants extracted from PDB files.
    """
    def __init__(self):
        self.pdb_file = None
        self.pdb = None
        self.molecule_type = None
        self.chain = None
        self.resnum = None
        self.resname = None
        self.content = None
        self.atoms = {}
        self.rdkit_mol = None
        self.pattern = None
        self.reactant_atoms = []

    def load_ligand(self, pdb_input, ligand_chain="L"):
        self.pdb = PDB(pdb_input)
        self.pdb_file = pdb_input
        self.chain = ligand_chain
        self.pdb.read_conect() # We need conectivity to get BONDS
        ligand = self.pdb.get_ligand(ligand_chain)
        self.resnum = int(get_resnum_from_line(ligand.split("\n")[0]))
        self.resname = get_resname_from_line(ligand.split("\n")[0])
        self.read_atoms(ligand)
        self.content = ligand + "\n" + "".join(self.pdb.conect_section) # List -> Str
        print("Ligand:")
        print(self.content)
        self.molecule_type = "ligand"

    def load_residue(self, pdb_input, residue_number, residue_chain):
        self.pdb = PDB(pdb_input)
        self.pdb.read_conect()
        self.pdb_file = pdb_input
        self.chain = residue_chain
        self.resnum = residue_number
        residue_content = self.pdb.get_residue(residue_number, residue_chain)
        self.read_atoms(residue_content)
        self.content = residue_content + "\n" + "".join(self.pdb.conect_section)
        self.resname = get_resname_from_line(residue_content.split("\n")[0])
        print("Residue:")
        print(self.content)
        self.molecule_type = "residue"

    def read_atoms(self, pdb_object):
        pdb_lines = pdb_object.split("\n")
        for n, line in enumerate(pdb_lines):
            pdb_at_name = get_atom_pdb_name_from_line(line).strip()
            self.atoms[pdb_at_name] = n

    def search_reacting_atoms(self):
        if not self.reactant_atoms: # To avoid multiplicity
            self.__transform_to_rdkit()
            if type(self.pattern) is str: 
                patt = Chem.MolFromSmarts(self.pattern)
                atom_indexes = self.rdkit_mol.GetSubstructMatch(patt)
                for at_id in atom_indexes:
                    self.reactant_atoms.append(get_specifyc_atom_pdb_name(self.content, 
                                                                          at_id))
            if type(self.pattern) is list:
                for pattern in self.pattern:
                    patt = Chem.MolFromSmarts(pattern)
                    atom_indexes = self.rdkit_mol.GetSubstructMatch(patt)
                    if atom_indexes:
                        for at_id in atom_indexes:
                            self.reactant_atoms.append(
                                                get_specifyc_atom_pdb_name(self.content, 
                                                                           at_id)
                                                      )
            if not self.reactant_atoms:
                raise ValueError("Reaction pattern not found in reactant")

    def get_connection_tree(self):
        self.__transform_to_rdkit()
        mol = self.rdkit_mol
        connections = []
        for bond in mol.GetBonds():
            connections.append([get_specifyc_atom_pdb_name(self.content, 
                                                           bond.GetBeginAtomIdx()),
                                get_specifyc_atom_pdb_name(self.content, 
                                                           bond.GetEndAtomIdx())])
        return connections

    def __transform_to_rdkit(self):
        if not self.content:
            raise AttributeError("The content is empty. Check that you have already load a ligand or a residue before transform it!")
        self.__create_temporary()
        self.rdkit_mol = Chem.MolFromPDBFile("tmp.pdb", removeHs=False)
        self.__delete_temporary()

    def __create_temporary(self):
        create_file(self.content, "tmp.pdb")

    def __delete_temporary(self):
        os.remove("tmp.pdb")
