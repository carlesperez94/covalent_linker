import os
from offpele.topology import Molecule
from offpele.template import Impact
from offpele.utils import get_data_file_path
import frag_pele.Covalent.correct_template_of_backbone_res as cov
from covalent_linker.constants import SCH_PATH, DATA

class ResidueLigandTemplate:

    def __init__(self, ligand_pdb, aminoacid, data_folder=DATA):
        self.ligand_pdb = ligand_pdb
        self.aminoacid = aminoacid
        self.data_folder = data_folder
        self.forcefield = 'OPLS2005'
        self.sch_path = SCH_PATH

    def set_OPLS2005_forcefield(self):
        self.forcefield = 'OPLS2005'

    def set_OFF_forcefield(self):
        self.forcefield = 'OpenFF'

    def create_aa_template_path(self):
        if self.forcefield == 'OPLS2005':
            path = os.path.join(self.data_folder, 
                                'Templates/OPLS2005/Protein',
                                self.aminoacid.lower())
        if self.forcefield == 'OpenFF':
            path = os.path.join(self.data_folder, 
                                'Templates/OPLS2005/Protein',
                                self.aminoacid.lower()) # Aminoacids are not parametrized in OFF yet.
        return path

    def get_template(self, outpath='grw'):
        os.environ['SCHRODINGER'] = self.sch_path
        m = Molecule(self.ligand_pdb, 
                     core_constraints=[' CA ', ' C  ', ' N  ']) 
        m.parameterize(self.forcefield)
        impact = Impact(m)
        impact.write(outpath)
        aa_template = self.create_aa_template_path()  
        cov.correct_template(outpath, aa_template)
        print("Template modified in {}.".format(outpath))

 
