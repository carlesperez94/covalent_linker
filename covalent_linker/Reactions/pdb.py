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
        self.read_all()

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
                print("Deleting... {}".format(line))
                self.lines.remove(line)
        self.__update_content_from_lines()

    def delete_residue(self, chain, resnum):
        for line in self.atom_section:
            chain_l = get_chain_from_line(line).strip()
            resnum_l = get_resnum_from_line(line).strip()
            if str(chain_l) == str(chain) and int(resnum_l) == int(resnum):
                print("Deleting... {}".format(line))
                try:
                    self.lines.remove(line)
                except ValueError:
                    print("{} does not exist".format(line))
        self.__update_content_from_lines()

    def replace_residue(self, chain, resnum,  residue_new):
        new_residue_lines = residue_new.split("\n")
        for n, line in enumerate(self.lines):
            if line.startswith("ATOM"):
                chain_l = get_chain_from_line(line).strip()
                resnum_l = get_resnum_from_line(line).strip()
                if str(chain_l) == str(chain) and int(resnum_l) == int(resnum):
                    self.delete_residue(chain, resnum)
                    for new_l in new_residue_lines:
                        self.lines.insert(n, new_l+"\n")
                        n = n+1
                    break
        self.__update_content_from_lines()

    def join_ligand_to_residue(self, res_chain, resnum, ligand_chain="L"):
        residue = self.get_residue(resnum, res_chain)
        ligand = self.get_ligand(ligand_chain)
        ligand_lines = ligand.split("\n")
        new_lig_lines = []
        for line in ligand_lines:
            if line.startswith("HETATM"):
                line = replace_hetatm_to_atom_in_line(line) 
                line = replace_chain_in_line(line, res_chain) 
                line = replace_resnum_in_line(line, resnum)
                new_lig_lines.append(line)
        ligname = get_resname_from_line(new_lig_lines[0])
        residue_lines = residue.split("\n")
        new_res_lines = []
        for line in residue_lines:
            if line.startswith("ATOM"):
                line = replace_resname_in_line(line, ligname)
                new_res_lines.append(line)
        new_res = "\n".join(new_res_lines) + "\n" + "\n".join(new_lig_lines)
        lignum = int(get_resnum_from_line(ligand_lines[0]))
        self.delete_residue(ligand_chain, lignum)
        self.replace_residue(res_chain, resnum, new_res)
    
    def __update_content_from_lines(self):
        self.content = "".join(self.lines)

    def write_content(self, outfile):
        with open(outfile, "w") as out:
            out.write(self.content)

def replace_hetatm_to_atom_in_line(line):
    line = list(line)
    line[0:6] = "ATOM  "
    line = "".join(line)
    return line

def replace_chain_in_line(line, new_chain):
    line = list(line)
    line[21] = new_chain
    line = "".join(line)
    return line

def replace_resnum_in_line(line, resnum):
    line = list(line)
    line[22:26] = "{:4d}".format(resnum)
    line = "".join(line)
    return line

def replace_resname_in_line(line, resname):
    line = list(line)
    line[17:20] = "{:3s}".format(resname)
    line = "".join(line)
    return line

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



