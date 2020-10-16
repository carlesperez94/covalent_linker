import os
from rdkit import Chem
from lib_prep.FragmentTools import tree_detector
from frag_pele.Growing import add_fragment_from_pdbs

from fragment_keys import fragments_relations
from pdb import PDB, get_specifyc_atom_pdb_name, get_resnum_from_line, get_atom_pdb_name_from_line


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
        self.read_atoms(self.pdb.get_residue(residue_number, residue_chain))
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

    def read_atoms(self, pdb_object):
        pdb_lines = pdb_object.split("\n")
        for n, line in enumerate(pdb_lines):
            pdb_at_name = get_atom_pdb_name_from_line(line).strip()
            self.atoms[pdb_at_name] = n

    def search_reacting_atoms(self):
        if not self.reactant_atoms: # To avoid multiplicity
            self.transform_to_rdkit()
            patt = Chem.MolFromSmarts(self.pattern) 
            atom_indexes = self.rdkit_mol.GetSubstructMatch(patt)
            for at_id in atom_indexes:
                self.reactant_atoms.append(get_specifyc_atom_pdb_name(self.content, at_id))

    def get_connection_tree(self):
        self.transform_to_rdkit()
        mol = self.rdkit_mol
        connections = []
        for bond in mol.GetBonds():
            connections.append([get_specifyc_atom_pdb_name(self.content, bond.GetBeginAtomIdx()),
                                get_specifyc_atom_pdb_name(self.content, bond.GetEndAtomIdx())])
        return connections     

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

    def get_bond_previous_bond_to_reaction(self):
        if not self.reaction_type:
            raise ValueError("Reaction Type: None. Please, select a reaction type!")
        connections = self.reactants[0].get_connection_tree()
        self.get_atoms_from_reaction()
        # Getting atoms from ligand
        self.reactants[0].search_reacting_atoms()
        reaction_bond = self.reactants[0].reactant_atoms
        # We first select the first (the most proximal to the ligand part) atom .
        first_atom_of_reaction_bond = reaction_bond[0]
        # Then. we must select which bond is connecting with this atom to get the previous bond to the reacting one.
        for bond in connections:
            if bond[1] == first_atom_of_reaction_bond:
                return bond[0].strip(), bond[1].strip()

    def prepare_intermediate(self):
        core_name, terminal_name = self.get_bond_previous_bond_to_reaction()
        fr_pdb, fr_core_name, fr_core_h = fragments_relations[self.reaction_type]
        outpath = os.path.join(get_directory_from_filepath(self.reactants[0].pdb_file), "reaction")
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        add_fragment_from_pdbs.main(self.reactants[0].pdb_file, 
                                    fr_pdb,
                                    pdb_atom_core_name=core_name,
                                    pdb_atom_fragment_name=fr_core_name, 
                                    steps=0,
                                    core_chain=self.reactants[0].chain,
                                    fragment_chain="L",
                                    output_file_to_tmpl="intermediate.pdb",
                                    output_file_to_grow="initial.pdb",
                                    h_core=terminal_name,
                                    h_frag=fr_core_h,  
                                    rename=False,
                                    threshold_clash=1.70, 
                                    output_path=outpath, 
                                    only_grow=False)
        intermediate_path = os.path.join(outpath, "pregrow", "intermediate.pdb")
        return intermediate_path

    def add_intermediate(self, intermediate_path):
        replaced_atom, static_atom = self.get_bond_previous_bond_to_reaction()
        outpath = os.path.join(get_directory_from_filepath(self.reactants[0].pdb_file), "reaction")
        # Add the intermediate (fragment) onto the COMPLEX PDB
        add_fragment_from_pdbs.main(self.reactants[0].pdb_file, 
                                    intermediate_path,
                                    pdb_atom_core_name=self.reactants[1].reactant_atoms[0].strip(), # Atom of the reactant 2
                                    pdb_atom_fragment_name=static_atom, 
                                    steps=0,
                                    core_chain=self.reactants[1].chain,
                                    fragment_chain=self.reactants[0].chain,
                                    output_file_to_tmpl="product.pdb",
                                    output_file_to_grow="initial.pdb",
                                    h_core=None,
                                    h_frag=None, 
                                    rename=False,
                                    threshold_clash=1.70, 
                                    output_path=outpath, 
                                    only_grow=False,
                                    cov_res="{}:{}".format(self.reactants[1].chain, self.reactants[1].resnum))
        print("{} completed successfully.".format(self.reaction_type))

    def apply_reaction(self):
         intermediate = self.prepare_intermediate()
         self.add_intermediate(intermediate)


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

def get_directory_from_filepath(path):
    return os.path.dirname(path)    
