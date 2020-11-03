import os

from rdkit import Chem
from lib_prep.FragmentTools import tree_detector
from frag_pele.Growing import add_fragment_from_pdbs
from frag_pele.Helpers.plop_rot_temp import prepare_pdb

from covalent_linker.Reactions.fragment_keys import fragments_relations
from covalent_linker.Reactions.pdb import PDB, get_specifyc_atom_pdb_name, get_resnum_from_line, get_atom_pdb_name_from_line
from covalent_linker.Reactions.reactants import Reactant
from covalent_linker.Helpers.filemanager import get_directory_from_filepath 
from covalent_linker.constants import SCH_PATH 

class Reaction:
    """
    Pattern class to define reactions.
    """
    def __init__(self):
        self.reactants = []
        self.reaction_type = None
        self.hydrogens_to_rm = None
        self.pdb_to_link = None

    def set_nucleophilic_addition_triple_bond_reaction(self):
        if len(self.reactants) != 2:
            raise ValueError("You need two reactants in this reaction")
        self.reaction_type = "Nucleophilic addition to triple bond"
        self.reactants[0].pattern = "[CX2]#[NX1]"
        self.reactants[1].pattern = ["[#16X2H]", "[OX2H]"] # Thiol, hydroxil

    def get_atoms_from_reaction(self):
        for react in self.reactants:
            react.search_reacting_atoms()
            print("{} reaction atoms for {} {}".format(self.reaction_type,
                                                       react.molecule_type, 
                                                       react.reactant_atoms))

    def apply_reaction(self, output_file="output.pdb", interm_out="intermediate.pdb" , free_out="product_free.pdb",
                       complex_out="product_complex.pdb",
                       outpath=os.path.join(os.getcwd(), "reaction")):
        intermediate = self.__prepare_intermediate(interm_out, complex_out, outpath)
        self.__add_intermediate(intermediate, free_out, "bad.pdb", outpath)
        self.__prepare_pdb_after_reaction()
        self.__write_pdb(outpath, output_file)
        self.__correct_pdb_output(output_file, outpath)


    def __prepare_pdb_after_reaction(self):
        if not self.hydrogens_to_rm:
            raise ValueError("Hydrogens to remove were not set. Please, apply a reaction before using this function")
        for num, hydrogenreact in enumerate(zip(self.hydrogens_to_rm, self.reactants)):
            hydrogen, react = hydrogenreact
            print("Deleting {}...".format(hydrogen.name))
            if num == 0:
                resnum = 1
            else:
                resnum = react.resnum
            self.pdb_to_link.delete_atom(chain=react.chain,
                                         resnum=resnum,
                                         atom_name=hydrogen.name)

    def __correct_pdb_output(self, output_file="output.pdb",
                             outpath=os.path.join(os.getcwd(), "reaction")):
        pdb = PDB(os.path.join(outpath, output_file))
        ch_r = self.reactants[1].chain
        ch_l = self.reactants[0].chain
        nu_r = self.reactants[1].resnum
        pdb.join_ligand_to_residue(ch_r, nu_r, ch_l)
        pdb.write_content(os.path.join(outpath, output_file))

    def __write_pdb(self, outpath, outfile):
        self.pdb_to_link.write_content(os.path.join(outpath, outfile))

    def __get_bond_previous_bond_to_reaction(self):
        if not self.reaction_type:
            raise ValueError("Reaction Type: None. Please, select a reaction type!")
        connections = self.reactants[0].get_connection_tree()
        self.get_atoms_from_reaction()
        # Getting atoms from ligand
        self.reactants[0].search_reacting_atoms()
        reaction_bond = self.reactants[0].reactant_atoms
        # We first select the first (the most proximal to the ligand part) atom .
        first_atom_of_reaction_bond = reaction_bond[0].strip()
        # Then. we must select which bond is connecting with this atom to get the previous bond to the reacting one.
        for bond in connections:
            if reaction_bond[0] in bond and reaction_bond[1] not in bond:
                prev_bond = [bond[0].strip(), bond[1].strip()]
                if prev_bond[0] == first_atom_of_reaction_bond: # Check if the first atom of the bond is or not the proximal 
                    prev_bond = prev_bond[::-1] # If is not, change the order to get the correct core-fragment direction
                return prev_bond

    def __prepare_intermediate(self, free_file="intermediate.pdb", complex_file="product_complex.pdb",
                               outpath=os.path.join(os.getcwd(), "reaction")):
        core_name, terminal_name = self.__get_bond_previous_bond_to_reaction()
        fr_pdb, fr_core_name, fr_core_h = fragments_relations[self.reaction_type]
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        add_fragment_from_pdbs.main(self.reactants[0].pdb_file, 
                                    fr_pdb,
                                    pdb_atom_core_name=core_name,
                                    pdb_atom_fragment_name=fr_core_name, 
                                    steps=0,
                                    core_chain=self.reactants[0].chain,
                                    fragment_chain="L",
                                    output_file_to_tmpl=free_file,
                                    output_file_to_grow=complex_file,
                                    h_core=terminal_name,
                                    h_frag=fr_core_h,  
                                    rename=False,
                                    threshold_clash=1.70, 
                                    output_path=outpath, 
                                    only_grow=False)
        # 'pregrow' folder is created by FragPELE
        self.pdb_to_link = PDB(os.path.join(outpath, "pregrow", complex_file))
        self.pdb_to_link.read_all()
        intermediate_path = os.path.join(outpath, "pregrow", free_file)
        return intermediate_path

    def __add_intermediate(self, intermediate_path, free_file="product_free.pdb", 
                           complex_file="bad.pdb",
                           outpath=os.path.join(os.getcwd(), "reaction")):
        replaced_atom, static_atom = self.__get_bond_previous_bond_to_reaction()
        # Add the intermediate (fragment) onto the COMPLEX PDB
        out_add = add_fragment_from_pdbs.main(self.reactants[0].pdb_file, 
                                              intermediate_path,
                                              pdb_atom_core_name=self.reactants[1].reactant_atoms[0].strip(), # Atom of the reactant 2
                                              pdb_atom_fragment_name=static_atom, 
                                              steps=0,
                                              core_chain=self.reactants[1].chain,
                                              fragment_chain=self.reactants[0].chain,
                                              output_file_to_tmpl=free_file,
                                              output_file_to_grow=complex_file,
                                              h_core=None,
                                              h_frag=None, 
                                              rename=False,
                                              threshold_clash=1.70, 
                                              output_path=outpath, 
                                              only_grow=False,
                                              cov_res="{}:{}".format(self.reactants[1].chain, 
                                                                     self.reactants[1].resnum))
        os.remove("bad.pdb") # Clear unnecessary files from add_fragments...
        os.remove("merged.pdb")
        self.hydrogens_to_rm = out_add[1][::-1] # Reversed to fit with the reactants
        curr_dir = os.getcwd()
        os.chdir(os.path.join(outpath, "pregrow"))
        prepare_pdb(pdb_in=free_file, 
                    pdb_out=free_file, 
                    sch_path=SCH_PATH)
        os.chdir(curr_dir)
        print("{} completed successfully.".format(self.reaction_type))

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
        


