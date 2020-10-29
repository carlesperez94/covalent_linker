import os
from offpele.topology import Molecule
from offpele.template import Impact
from offpele.utils import get_data_file_path
import frag_pele.Covalent.correct_template_of_backbone_res as cov

class ResidueLigandTemplate:

    def __init__(self, ligand_pdb):
        self.ligand_pdb = ligand_pdb
        self.forcefield = 'OPLS2005'
        self.sch_path = '/gpfs/projects/bsc72/SCHRODINGER_ACADEMIC'

    def set_OPLS2005_forcefield():
        self.forcefield = 'OPLS2005'

    def set_OFF_forcefield():
        self.forcefield = 'OpenFF'

    def get_template(self, outpath='grw'):
        os.environ['SCHRODINGER'] = self.sch_path
        m = Molecule(self.ligand_pdb, 
                     core_constraints=[' CA ', ' C  ', ' N  ']) 
        m.parameterize(self.forcefield)
        impact = Impact(m)
        impact.write(outpath)  
        cov.correct_template(outpath) 

 
