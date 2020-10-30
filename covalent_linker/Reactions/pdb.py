class PDB:
    def __init__(self, pdb_file):
        """
        Class to read and process PDB files with a single model.
        :param in_pdb: pdb file path
        :type pdb_file: str
        """
        self.pdb_file = pdb_file
        self.content = None
        self.lines = None
        self.atom_section = None
        self.conect_section = None

    def read_content(self):
        """
        Reads the content of the PDB file.
        """
        with open(self.pdb_file) as infile:
            self.content = infile.read()

    def read_lines(self):
        """
        Reads the content of the PDB file.
        """
        with open(self.pdb_file) as infile:
            self.lines = infile.readlines()

    def read_atoms_section(self):
        """
        Reads the ATOM section of the PDB as a list.
        """
        self.atom_section = []
        if not self.lines:
            self.read_lines()
        for line in self.lines:
            if "ATOM" in line[0:4] or "HETATM" in line[0:6]:
                self.atom_section.append(line)

    def read_conect(self):
        """
        Reads the CONECT section of a PDB as a list.
        """
        self.conect_section = []
        if not self.lines:
            self.read_lines()
        for line in self.lines:
            if "CONECT" in line[0:6]:
                self.conect_section.append(line)

    def get_chain(self, chain):
        """
        Gets the lines (as str) of the selected chain.
        """
        if not self.atom_section:
            self.read_atoms_section()
        chain_lines = []
        for at_line in self.atom_section:
            if at_line[21:22] == chain:
                chain_lines.append(at_line)
        return "".join(chain_lines)

    def get_ligand(self, ligand_chain="L"):
        """
        Gets the ligand chain and checks that all lines are heteroatoms.
        """
        ligand = self.get_chain(ligand_chain)
        for line in ligand.split("\n"):
            if not "HETATM" in line[0:6] and not line == "":
                print("Selected chain: {}".format(ligand_chain))
                print("Line:")
                print(line)
                raise TypeError("The selected chain does not contain HETATM in all lines!")
        return ligand

    def get_residue(self, resnum, chain):
        """
        Gets the lines (as str) of the specified residue number and chain.
        """
        res_chain = self.get_chain(chain)
        residue = []
        for line in res_chain.split("\n"):
            if str(resnum) == str(line[22:26].strip()):
                residue.append(line)
        return "\n".join(residue)

    def read_all(self):
        self.read_content()
        self.read_lines()
        self.read_atoms_section()
        self.read_conect()
        

    def delete_atom(self, chain, resnum, atom_name):
        for line in self.atom_section:
            chain_l = get_chain_from_line(line).strip()
            resnum_l = get_resnum_from_line(line).strip()
            name_l = get_atom_pdb_name_from_line(line).strip()
            if str(chain_l) == str(chain) and int(resnum_l) == int(resnum) and str(name_l) == str(atom_name):
                print(line)
                self.lines.remove(line)
    
    def update_content_from_lines(self):
        self.content = "".join(self.lines)

    def write_content(self, outfile):
        with open(outfile, "w") as out:
            out.write(self.content)


def get_chain_from_line(line):
    return line[21]

def get_resnum_from_line(line):
    return line[22:26]

def get_atom_pdb_name_from_line(line):
    return line[12:16]

def get_resname_from_line(line):
    return line[17:20]

def get_specifyc_atom_pdb_name(content, line_idx):
    lines = content.split("\n") 
    return get_atom_pdb_name_from_line(lines[line_idx])



