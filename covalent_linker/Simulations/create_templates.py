import os
from peleffy.topology import Molecule, RotamerLibrary
from peleffy.template import Impact
from peleffy.utils import get_data_file_path
import frag_pele.Covalent.correct_template_of_backbone_res as cov
from frag_pele.Helpers import folder_handler
from covalent_linker.constants import SCH_PATH, DATA

class ResidueLigandTemplate:

    def __init__(self, ligand_pdb, aminoacid, data_folder=DATA):
        self.ligand_pdb = ligand_pdb
        self.aminoacid = aminoacid
        self.data_folder = data_folder
        self.__forcefield = 'OPLS2005'
        self.sch_path = SCH_PATH

    def set_OPLS2005_forcefield(self):
        self.__forcefield = 'OPLS2005'

    def set_OFF_forcefield(self):
        self.__forcefield = 'OpenFF'

    def __create_aa_template_path(self):
        if self.__forcefield == 'OPLS2005':
            path = os.path.join(self.data_folder, 
                                'Templates/OPLS2005/Protein',
                                self.aminoacid.lower())
        if self.__forcefield == 'OpenFF':
            path = os.path.join(self.data_folder, 
                                'Templates/OPLS2005/Protein',
                                self.aminoacid.lower()) # Aminoacids are not parametrized in OFF yet.
        return path

    def __get_template_and_rot(self, template_path='grw', rot_path='GRW.rot.assign'):
        os.environ['SCHRODINGER'] = self.sch_path
        m = Molecule(self.ligand_pdb, 
                     core_constraints=[' CA ', ' C  ', ' N  ']) 
        m.parameterize(self.__forcefield)
        impact = Impact(m)
        impact.write(template_path)
        aa_template = self.__create_aa_template_path()  
        cov.correct_template(template_path, aa_template)
        print("Template modified in {}.".format(template_path))
        rotamer_library = RotamerLibrary(m)
        rotamer_library.to_file(rot_path)
        print("Rotamer library stored in {}".format(rot_path))

    def get_datalocal(self, outdir=".", name="grw"):
        folder_handler.check_and_create_DataLocal(working_dir=outdir)
        if self.__forcefield == 'OPLS2005':
            datalocal_temp_path = os.path.join(outdir,  
                                          "DataLocal/Templates/OPLS2005/Protein",
                                          name)
        if self.__forcefield == 'OpenFF':
            datalocal_temp_path = os.path.join(outdir,  
                                          "DataLocal/Templates/OPLS2005/Protein",
                                          name)
        datalocal_rot_path = os.path.join(outdir, "DataLocal/LigandRotamerLibs", 
                                          "{}.rot.assign".format(name.upper()))
        self.__get_template_and_rot(template_path=datalocal_temp_path, 
                                    rot_path=datalocal_rot_path)

 
