import os
import shutil
import subprocess
from frag_pele.Helpers import constraints
import covalent_linker.Helpers.template_builder as tb
from covalent_linker.constants import LICENSE

class PELEControlFile:

    def __init__(self, template, outpath="pele.conf", license=LICENSE, 
                 results_path="simulation_output", steps=100, overlap=0.7, 
                 chain="L", con_backbone=5, con_ter=5, con_interval=10,
                 center=[0,0,0], temperature=1000, seed=1000, steering=0, 
                 translation_high=0.2, translation_low=0.1, rotation_high=0.1, 
                 rotation_low=0.05, radius=6, reschain=None, resnum=None):
        self.template = template
        self.pdbs = []
        self.license = license
        self.results_path = results_path
        self.steps = steps
        self.overlap = overlap
        self.chain = chain
        self.center = center
        self.temperature = temperature
        self.seed = seed
        self.steering = steering
        self.translation_high = translation_high
        self.translation_low = translation_low
        self.rotation_high = rotation_high
        self.rotation_low = rotation_low
        self.con_backbone = con_backbone
        self.con_ter = con_ter
        self.con_interval = con_interval
        self.radius = radius
        self.reschain = reschain
        self.resnum = resnum
        self.pdb_lines = None
        self.constraints = None
        self.content = None
 

    def load_pdbs(self, pdb_files):
        if type(pdb_files) == str:
            self.pdbs.append(pdb_files)
        if type(pdb_files) == list:
            for pdb in pdb_files:
                self.pdbs.append(pdb)
        self.__create_pdb_lines()
        self.constraints = "\n".join(constraints.retrieve_constraints(self.pdbs[0], 
                                                                      {}, {},
                                                                      self.con_backbone, 
                                                                      self.con_ter, 
                                                                      self.con_interval))
   
    def __create_pdb_lines(self):
        self.pdb_lines = None
        lines = []
        for pdb in self.pdbs:
            line = '{"files" : [{"path": "%s" }] }' % (pdb)
            lines.append(line)
        self.pdb_lines = ",\n".join(lines)
     
    def fill_content(self):
        # Definition of the keywords that we are going to substitute from the template
        keywords = {"LICENSE": self.license,
                    "RESULTS_PATH": self.results_path,
                    "CHAIN": self.chain,
                    "CONSTRAINTS": self.constraints,
                    "CENTER": self.center, 
                    "PDB": self.pdb_lines,
                    "STEPS": self.steps,
                    "OVERLAP": self.overlap,
                    "TEMPERATURE": self.temperature,
                    "SEED": self.seed,
                    "STEERING": self.steering,
                    "TRANSLATION_HIGH": self.translation_high,
                    "TRANSLATION_LOW": self.translation_low,
                    "ROTATION_HIGH": self.rotation_high,
                    "ROTATION_LOW": self.rotation_low,
                    "RADIUS": self.radius,
                    "RESCHAIN": self.reschain,
                    "RESNUM": self.resnum
                    }
        template_filled = tb.TemplateBuilder(self.template, keywords)
        self.content = template_filled.conf



        
