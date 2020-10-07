import os
from pdb import PDB
from rdkit import Chem


class Reactant:
    """
    It defines reactants extracted from PDB files.
    """
    def __init__(self):
        self.pdb_file = None
        self.pdb = None
        self.molecule_type = None
        self.content = None
        self.rdkit_mol = None
    
    def load_ligand(self, pdb_input, ligand_chain="L"):
        self.pdb = PDB(pdb_input)
        self.pdb_file = pdb_input
        self.content = self.pdb.get_ligand(ligand_chain)
        print("Ligand:")
        print(self.content)
        self.molecule_type = "ligand"

    def load_residue(self, pdb_input, residue_number, residue_chain):
        self.pdb = PDB(pdb_input)
        self.pdb_file = pdb_input
        self.content = self.pdb.get_residue(residue_number, residue_chain)
        print("Residue:")
        print(self.content)
        self.molecule_type = "residue"

    def transform_to_rdkit(self):
        if not self.content:
            raise AttributeError("The content is empty. Check that you have already load a ligand or a residue before transform it!")
        tmp_filename = "tmp.pdb"
        create_temporary_file(self.content, tmp_filename)
        self.rdkit_mol = Chem.MolFromPDBFile(tmp_filename, removeHs=False)
        os.remove(tmp_filename)


class Reaction:
    """
    Pattern class to define reactions.
    """
    def __init__(self):
        self.reactant1 = None
        self.reactant2 = None
        self.reaction_type = None

class LigandResidueReaction(Reaction):
    """
    Class to define reactions between ligands and residues (covalent linking).
    """
    def __init__(self):
        super().__init__()

    def load_pdb(self, pdb_file, residue_chain, residue_num, ligand_chain="L"):
        self.reactant1 = Reactant()
        self.reactant1.load_ligand(pdb_file, ligand_chain)
        self.reactant2 = Reactant()
        self.reactant2.load_residue(pdb_file, residue_num, residue_chain) 

def create_temporary_file(content, filename):
    with open(filename, "w") as out:
        out.write(content)


