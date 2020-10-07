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
        self.pattern = None
        self.reactant_atoms = None
    
    def load_ligand(self, pdb_input, ligand_chain="L"):
        self.pdb = PDB(pdb_input)
        self.pdb_file = pdb_input
        self.pdb.read_conect() # We need conectivity to get BONDS
        self.content = self.pdb.get_ligand(ligand_chain) + "\n" + "".join(self.pdb.conect_section) # List -> Str
        print("Ligand:")
        print(self.content)
        self.molecule_type = "ligand"

    def load_residue(self, pdb_input, residue_number, residue_chain):
        self.pdb = PDB(pdb_input)
        self.pdb.read_conect()
        self.pdb_file = pdb_input
        self.content = self.pdb.get_residue(residue_number, residue_chain) + "\n" + "".join(self.pdb.conect_section)
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

    def search_pattern(self):
        self.transform_to_rdkit()
        patt = Chem.MolFromSmarts(self.pattern) 
        self.reactant_atoms = self.rdkit_mol.GetSubstructMatch(patt)


class Reaction:
    """
    Pattern class to define reactions.
    """
    def __init__(self):
        self.reactants = []
        self.reaction_type = None

    def set_triple_bond_addition_reaction(self):
        if len(self.reactants) != 2:
            raise ValueError("You need two reactants in this reaction")
        self.reaction_type = "Triple bond addition"
        self.reactants[0].pattern = "[NX1]#[CX2]"
        self.reactants[1].pattern = "[#16X2H]"

    def get_atoms_from_reaction(self):
        for react in self.reactants:
            react.search_pattern()
            print("Reaction atoms for {} {}".format(react.molecule_type, 
                                                    react.reactant_atoms))


class LigandResidueReaction(Reaction):
    """
    Class to define reactions between ligands and residues (covalent linking).
    """
    def __init__(self):
        super().__init__()

    def load_pdb(self, pdb_file, residue_chain, residue_num, ligand_chain="L"):
        r1 = Reactant()
        r1.load_ligand(pdb_file, ligand_chain)
        self.reactants.append(r1)
        r2 = Reactant()
        r2.load_residue(pdb_file, residue_num, residue_chain) 
        self.reactants.append(r2)
        

def create_temporary_file(content, filename):
    with open(filename, "w") as out:
        out.write(content)




