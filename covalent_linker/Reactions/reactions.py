from rdkit import Chem

class Reaction:
    def __init__(self):
        self.pdb_file = None
        self.ligand_chain = None
        self.ligand = None
        self.residue_id = None
        self.residue = None

    def get_ligand(self):
        return self.ligand

    def get_residue(self):
        return self.residue

