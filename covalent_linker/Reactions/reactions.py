import os
from rdkit import Chem
from lib_prep.FragmentTools import tree_detector
from frag_pele.Growing import add_fragment_from_pdbs

from fragment_keys import fragments_relations
from pdb import PDB, get_specifyc_atom_pdb_name, get_resnum_from_line


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
        self.content = None
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
        self.content = self.pdb.get_residue(residue_number, residue_chain) + "\n" + "".join(self.pdb.conect_section)
        print("Residue:")
        print(self.content)
        self.molecule_type = "residue"

    def transform_to_rdkit(self):
        if not self.content:
            raise AttributeError("The content is empty. Check that you have already load a ligand or a residue before transform it!")
        self.create_temporary()
        self.rdkit_mol = Chem.MolFromPDBFile("tmp.pdb", removeHs=False)
        self.delete_temporary()

    def create_temporary(self):
        create_temporary_file(self.content, "tmp.pdb")

    def delete_temporary(self):
        os.remove("tmp.pdb")

    def search_reacting_atoms(self):
        self.transform_to_rdkit()
        patt = Chem.MolFromSmarts(self.pattern) 
        atom_indexes = self.rdkit_mol.GetSubstructMatch(patt)
        for at_id in atom_indexes:
            self.reactant_atoms.append(get_specifyc_atom_pdb_name(self.content, at_id))
         

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
        self.reactants[0].pattern = "[CX2]#[NX1]"
        self.reactants[1].pattern = "[#16X2H]"

    def get_atoms_from_reaction(self):
        for react in self.reactants:
            react.search_reacting_atoms()
            print("{} reaction atoms for {} {}".format(self.reaction_type,
                                                       react.molecule_type, 
                                                       react.reactant_atoms))

    def select_fragment(self):
        return fragments_relations[self.reaction_type]

    def apply_reaction(self):
        self.get_atoms_from_reaction()
        fr_pdb, fr_core_name, fr_core_h = self.select_fragment()
        print(self.reactants[0].pdb_file, fr_pdb, fr_core_name, fr_core_h)
        add_fragment_from_pdbs.main(self.reactants[0].pdb_file, fr_pdb, 
                                    pdb_atom_core_name=self.reactants[0].reactant_atoms[0], 
                                    pdb_atom_fragment_name=fr_core_name, steps=0, 
                                    core_chain=self.reactants[0].chain,
                                    fragment_chain="L", output_file_to_tmpl="growing_result.pdb", 
                                    output_file_to_grow="initialization_grow.pdb",
                                    h_core=self.reactants[0].reactant_atoms[1], 
                                    h_frag=fr_core_h, rename=False, 
                                    threshold_clash=1.70, output_path=".", only_grow=False)


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




