import glob
import shutil
import filecmp

from covalent_linker.Reactions.reactions import LigandResidueReaction
from covalent_linker.Reactions.pdb import PDB
from covalent_linker.Simulations.create_templates import ResidueLigandTemplate 

def test_nucleophilic_triple_CYS():
    r = LigandResidueReaction()
    r.load_pdb("input_files/test_cys_n_t.pdb", residue_chain="A", residue_num=145)
    r.set_nucleophilic_addition_triple_bond_reaction()
    r.apply_reaction(outpath="output_files/reaction", output_file="output.pdb")
    out1 = "output_files/reaction/output.pdb"
    exp1 = "expected_files/nucleoph_t_bond/output_cys.pdb"
    out2 = "output_files/reaction/pregrow/product_free.pdb"
    exp2 = "expected_files/nucleoph_t_bond/product_free_cys.pdb"
    assert filecmp.cmp(out1, exp1)
    assert filecmp.cmp(out2, exp2)

def test_template_creation():
    out1 = "output_files/DataLocal/Templates/OPLS2005/Protein/grw"
    exp1 = "expected_files/nucleoph_t_bond/grw"
    out2 = "output_files/DataLocal/LigandRotamerLibs/GRW.rot.assign"
    exp2 = "expected_files/nucleoph_t_bond/GRW.rot.assign"
    pdbfile = "output_files/reaction/pregrow/product_free.pdb"
    templ = ResidueLigandTemplate(pdbfile, "CYS")
    templ.get_datalocal(outdir="output_files", name="grw") 
    assert filecmp.cmp(out1, exp1)
    assert filecmp.cmp(out2, exp2)
    clear()

def test_nucleophilic_triple_SER():
    r = LigandResidueReaction()
    r.load_pdb("input_files/test_ser_n_t.pdb", residue_chain="A", residue_num=145)
    r.set_nucleophilic_addition_triple_bond_reaction()
    r.apply_reaction(outpath="output_files/reaction", output_file="output.pdb")
    out1 = "output_files/reaction/output.pdb"
    exp1 = "expected_files/nucleoph_t_bond/output_ser.pdb"
    out2 = "output_files/reaction/pregrow/product_free.pdb"
    exp2 = "expected_files/nucleoph_t_bond/product_free_ser.pdb"
    assert filecmp.cmp(out1, exp1)
    assert filecmp.cmp(out2, exp2)
    clear()

def clear():
    to_rm = glob.glob("output_files/*")
    for f in to_rm:
        shutil.rmtree(f)
