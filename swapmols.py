#!/usr/bin/env python3
'''
Script to modify pdb files --> Add Imidazolium to Cholesterol

Source molecules is the molecule to change
Target molecule is  the molecule to put into the structure

Example input:
    python3 morph.py -f bilayer.pdb -o output.pdb -s MOLNAME -t target_mo.pdb -i edit_specs.inp
'''

import argparse
import re
import numpy as np
import logging

LOGGER = logging.getLogger()
LOGGER.setLevel(logging.DEBUG)
FORMATTER = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
CH = logging.StreamHandler()
CH.setLevel(logging.DEBUG)
CH.setFormatter(FORMATTER)
LOGGER.addHandler(CH)

#### GLOBALS ####
PDB_REC = 'ATOM'
#                    Record   Serial   atomname    resn    resid        x       y           z    Rest
pdb_string_format = '{: <4}  {: >5} {: >2}{: <2}{}{: >4}{}{: >4}{}   {: >8.3f}{: >8.3f}{: >8.3f}{}'#{: >6.2f}{: >6.2f}'
pdb_pattern = r"""
        ^
        ATOM              (?# Record Type)
\s+
        (\d+)             (?# Serial numbers  Grp 1)
\s*
        (\w+?)([\d,\w,\']*)       (?# Atom names;     Grp 2+3)
\s+
        ([\w,\d]+)        (?# Residue name    Grp 4)
\s+
        (\d+)             (?# Resid           Grp 5)
\s+
        (-?\d+\.\d+)      (?# X               Grp 6-8)
\s*
        (-?\d+\.\d+)      (?# Y)
\s*
        (-?\d+\.\d+)      (?# Z)
        (.*)              (?# Remainder)
            """
pdb_pattern = ''.join(pdb_pattern.split()) # Remove all whitespaces in above defined pattern
regex_pdb = re.compile(pdb_pattern)

#### Parser arguments ####
parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store', nargs=1, metavar='input.pdb', required=True,
                    type=str, help='Input membrane structure file (pdb).')
parser.add_argument('-o', action='store', nargs=1, metavar='output.pdb', required=True,
                    type=str, help='Output membrane structure file (pdb).')
parser.add_argument('-t', action='store', nargs=1, metavar='structure_target.pdb', required=True,
                    type=str, help='Target molecule structure file (pdb).')
parser.add_argument('-i', action='store', nargs=1, metavar='edit_specifications.inp', required=True,
                    type=str, help='File that contains a list of all atoms to be edited.')
parser.add_argument('-s', action='store', nargs=1, metavar='<molname to be modifies>', required=True,
                    type=str, help='Name of molecule to change in source structure.\n'
                    'Syntax is described in function "read_input_atom_specs".')
parser.add_argument('-c', action='store', nargs=1, metavar='<conc/%>', required=False,
                    type=int, help='Concentration of changed molecule per bilayer')
args=parser.parse_args()   ### Arguments can be taken via args.<flagname>

INP_FNAME = args.f[0]
OUT_FNAME = args.o[0]
TARGSTRUCT_FNAME = args.t[0]
EDITSPECS_FNAME = args.i[0]
MOLNAME = args.s[0]
ATMSEXPL = []
try:
    CONCENTRATION = args.c[0]
except TypeError:
    CONCENTRATION = 100
####

LOGGER.info("Following input was given %s %s %s %s %s %s ",
            INP_FNAME, OUT_FNAME, TARGSTRUCT_FNAME, EDITSPECS_FNAME,
            MOLNAME, CONCENTRATION)

class transform_coordinates():
    ''' Class for handling all coordinate transformation calculations
        functions:
            -add_molecule: Add all atoms from target structure to source structure
            using
                -gen_transmat: Calculate transformation matrix from target to source system
                -calc_coords: Calculates new coordinates for target atom names
            -get_target_top_information: Gets molecule name, atomname order, and dictionary
                                         atom name to coordinate from target structure

    '''
    def __init__(self, refatmlist_src, refatmlist_targ,
                 skip_tar=None, mapping=None):
        '''
        refatm_l: List containing atomnames in order
                  specified as input arguments
        refatm_c: Dictionary containing atomname to coordinates of
                  atoms spanning internal coordinate system
        keep_src (optional): List of atomnames from target molecule
                             not to be modified
        skip_tar (optional): Linked to keep_src - List of atomnames that
                             are not to be added from target structure
        '''
        self.refatm_l_src  = refatmlist_src
        self.refatm_l_targ = refatmlist_targ

        if skip_tar is not None:
            self.skip_tar = skip_tar
        else:
            self.skip_tar = []
        if mapping is not None:
            self.mapping = mapping
        else:
            self.mapping = {}

        #
        returnvals = self.get_target_top_information(TARGSTRUCT_FNAME)
        self.atom_entry_order = returnvals[0]
        self.target_resname   = returnvals[1]
        self.atomname2relcoords_tar = returnvals[2]
        #
        self.number_of_atoms  = len(self.atom_entry_order)

    def add_molecule(self, fileobj, resid, atom_dict, other=''):
        ''' Main function to write all pdb lines for molecule <resid> in atom_entries variable
            refatom_dict: Dictionary for reference atoms forming coordinate system linked to their coordinates
            example: {'C1':np.array([xyz])}
        '''
        #print("Adding", resid)
        #print("Using", atom_dict.keys())
        transformation_matrix, origin_src = self.gen_transmat(atom_dict)
        atmname2coords = self.calc_coords(transformation_matrix, origin_src)
        for i in range(self.number_of_atoms):
            atmname = self.atom_entry_order[i]
            if atmname in self.skip_tar:
                xyz = atom_dict[self.mapping[atmname]]
            xyz = atmname2coords[atmname]
            xyz_str = '{: >8.3f}{: >8.3f}{: >8.3f}'.format(*xyz)
            pdbline = '{: <5} {: >5} {: >4}{}{: >4}{}{: >4}{}   {}{}'\
                .format(PDB_REC, i, atmname, ' ',
                self.target_resname.upper(), '', resid, ' ',
                xyz_str, other)
            print(pdbline, file=fileobj)

    def calc_coords(self, transmatrix, origin_src):
        ''' Calculates coordinates of atoms in target system
            using a transformationmatrix
        '''
        atmname2coords = {}
        for atm in self.atom_entry_order:
            xyz = self.atomname2relcoords_tar[atm]
            xyz_rel = np.array(np.dot(transmatrix, xyz))[0] + origin_src
            atmname2coords[atm] = xyz_rel
        return atmname2coords

    def gen_transmat(self, basis_coords_src):
        ''' Generates transformation matrix from target to source coordinate system
            The input parameters contain dictionaries with atmnames and coordinates
            of atoms forming coordinate system
            basis_coords_src: {atomname:<coords in source structure file>}
        '''
        LOGGER.debug("REFATM_L_SRC %s", self.refatm_l_src)
        e_src = []
        origin_src = basis_coords_src[self.refatm_l_src[0]]
        for atm in self.refatm_l_src[1:]:
            xyz = basis_coords_src[atm]
            vec = xyz - origin_src
            vec = vec/np.linalg.norm(vec)
            e_src.append(vec)
        LOGGER.debug("Matrix of e_src: \n %s", e_src)
        transmat = np.linalg.inv(np.matrix(e_src))
        LOGGER.debug("Translation matrix %s \n Origin coords %s", transmat, origin_src)
        return transmat, origin_src

    def get_target_top_information(self, topology_target):
        ''' Gets the order of atoms in a pdb structure file
            and returns:
                outlist        --> A list containing atomname entries
                                   in input file order
                resname        --> name of molecule
                atom2relcoords --> dictionary containing map of pdb
                                   atom name to relative coordinates
        '''
        atom2coords = {}
        refatms = {}
        outlist = []
        atom2relcoords = {}
        with open(topology_target, "r") as topf:
            for line in topf:
                regmatch = regex_pdb.match(line)
                if not regmatch:
                    continue
                grps = regmatch.groups() # [serial, atmtype, atmnumber, resname]
                atom1, atom2, resname = grps[1:4]
                atom = atom1+atom2
                xyz = np.array([float(x) for x in grps[5:8]])

                if atom in self.refatm_l_targ:
                    #print("Im here", atom)
                    refatms[atom] = xyz
                atom2coords[atom] = xyz
                outlist.append(atom)
        e_tar = []
        origin_tar = refatms[self.refatm_l_targ[0]]
        for atm in self.refatm_l_targ[1:]:
            xyz = refatms[atm]
            vec = xyz - origin_tar
            vec = vec/np.linalg.norm(vec)
            e_tar.append(vec)
        transmat = np.matrix(e_tar)
        for atm in outlist:
            xyz = atom2coords[atm]-origin_tar
            xyz_rel = np.array(np.dot(transmat, xyz))[0]
            atom2relcoords[atm] = xyz_rel
        return outlist, resname, atom2relcoords


def read_input_atom_specs(inputfilename):
    ''' Reads the input file specified in -i flag
        Requires 2 entries:
            source X1 X2 X3 X4  # 4 atoms that form coordinate system in source molecule
            target X1 X2 X3 X4  # Same as above for target molecule
        2 are optional:
            skip_tar    (see documentation of transform_coordinates class)
            keep_source
    '''
    outdict = {}
    with open(inputfilename, "r") as inpf:
        for line in inpf:
            nocomments = line.split('#')[0].strip()
            cols = nocomments.split()
            if not cols:
                continue
            specification = cols[0]
            pars = cols[1:]
            outdict[specification] = pars
    return outdict

def count_mols_of(molname, struct_file_source):
    ''' Count how many molecules of type <molname> are in source structure (specified in -f) '''
    reached_once = False
    counter = 0
    counter_atoms =0
    resid_old = None
    with open(struct_file_source, "r") as topf:
        for line in topf:
            regmatch = regex_pdb.match(line)
            if not regmatch:
                continue
            grps = regmatch.groups() # [serial, atmtype, atmnumber, resname]
            resname = grps[3]
            resid = grps[4]
            if resid_old is None:
                resid_old = resid
                continue
            if resname == molname:
                counter_atoms += 1
                if resid_old != resid:
                    reached_once = True
                    counter += 1
            resid_old = resid
    if not reached_once:
        counter += 1
    try:
        counter_atoms = int(counter_atoms/counter)
    except ZeroDivisionError:
        print("ERROR:",
              "The input molecule name seems not to be correct.",
              "No atoms found in source structure with this molecule name",
              )
        raise
    LOGGER.info("Mols found %s, atoms found: %s", counter, counter_atoms)
    return counter, counter_atoms

def check_input():
    ''' Checks if input is valid. Raise Error if not.'''
    if INP_FNAME == OUT_FNAME:
        raise ValueError("Input file must not be the same as output file:", INP_FNAME)
    elif TARGSTRUCT_FNAME == OUT_FNAME:
        raise ValueError("Input file must not be the same as output file:", TARGSTRUCT_FNAME)
    if INP_FNAME[-4:] != '.pdb' or TARGSTRUCT_FNAME[-4:] != '.pdb':
        print("WARNING are input structure files really in .pdb format?")



def main():
    '''
        Reads input pdb file, gather data of molecules to morph (in mol_pars) and passes data to class transform_coordinates
    '''
    check_input()
    inp_specs = read_input_atom_specs(EDITSPECS_FNAME)
    name_mapping = {}
    changed_mol_counter = 0
    unchanged_mol_counter = 0
    added_mols_total = 0
    count_atoms = 0
    n_mols, n_atoms_mol = count_mols_of(MOLNAME, INP_FNAME)
    n_mols_per_bilayer = n_mols//2
    if not n_mols_per_bilayer: # In the case that only one molecule is to be changed, counter, counter_atoms
        n_mols_per_bilayer = 1
    mols_to_change = int(n_mols_per_bilayer*CONCENTRATION*0.01)
    print("Morphing {} molecules per bilayer of {} molecules to get a concentration of {}%".format(mols_to_change, n_mols, CONCENTRATION))

    try:
        keep_src = inp_specs['keep_src']
        skip_tar = inp_specs['skip_tar']
        for k_src, s_tar in zip(keep_src, skip_tar):
            name_mapping[s_tar] = k_src
    except IndexError:
        keep_src = []
        skip_tar = []
    refatm_src = inp_specs['source']
    refatm_targ = inp_specs['target']

    struct_calc = transform_coordinates(refatm_src, refatm_targ, skip_tar=skip_tar, mapping=name_mapping)

    print('Reference target:{}\nReference source:{}\n'
           './{}\ntaking information from ./{}.\n\n'\
           .format(refatm_targ, refatm_src,
               INP_FNAME, TARGSTRUCT_FNAME))

    #### Start reading input and writing output bilayer structure ####
    with open(INP_FNAME, "r") as inpf, open(OUT_FNAME, "w") as outf:
        resid_old = 0     # To check if new residue sequence has begun
        resname_old = None
        first = True
        mol_pars = {}   # {<atmname>:XYZ}
        for line in inpf:

            regmatch = regex_pdb.match(line)
            if regmatch:
                grps = regmatch.groups()
                serial, atom_str, atom_int, resname, resid, x, y, z, other = grps
                atom = atom_str+atom_int    # letter and number of atomnames are handled differently

                # Start working on molecule or first molecule is already target:
                if resname_old == MOLNAME or (resname_old is None and resname == MOLNAME):
                    resid = int(resid)
                    xyz = np.array([float(val) for val in [x, y, z]])
                    mol_pars[atom] = xyz

                    # If enough molecules were morphed, but the residuals still have to be includeed:
                    if changed_mol_counter == mols_to_change and\
                            not (changed_mol_counter + unchanged_mol_counter == n_mols_per_bilayer):
                        if (count_atoms+1) / n_atoms_mol == 1:
                            LOGGER.info("Add source mol. Resid %s, Changed %s, Unchanged %s",
                                resid_old, changed_mol_counter, unchanged_mol_counter)
                            unchanged_mol_counter += 1
                            count_atoms = 0
                            if resname == MOLNAME:
                                print(line, file=outf, end='')
                        else:
                            print(line, file=outf, end='')
                            count_atoms += 1

                    # If all molecules were added
                    if changed_mol_counter + unchanged_mol_counter == n_mols_per_bilayer:
                        added_mols_total = added_mols_total + changed_mol_counter + unchanged_mol_counter
                        if unchanged_mol_counter != 0:# If not all residues are morphed
                            resid_old += 1
                        LOGGER.info("Reset. Nchanged %s and Nunchanged %s",
                            changed_mol_counter, unchanged_mol_counter)
                        changed_mol_counter = 0
                        unchanged_mol_counter = 0
                        LOGGER.info("N total added: %s", added_mols_total)
                        #continue

                    # If new resid entry begun --> pass output.pdb-file object
                    # First var prevents script from immediately adding mol if struct file begins with mol
                    if resid_old != resid and not first:
                        if changed_mol_counter < mols_to_change and added_mols_total < n_mols:
                            print("Add", resid_old, changed_mol_counter, unchanged_mol_counter)
                            struct_calc.add_molecule(outf, resid_old, mol_pars)
                            changed_mol_counter += 1
                            mol_pars = {}
                    if resid_old != resid and not first and resname != MOLNAME:
                        print(line, file=outf, end='')

                    first = False
                    resname_old = resname

                else:
                    print(line, file=outf, end='')

                resid_old = resid
            else:
                 #Last calculation and no other residues come after
                if resname_old == MOLNAME and added_mols_total < n_mols:
                    if changed_mol_counter < mols_to_change:
                        print("Added mols and totals", added_mols_total, n_mols)
                        print("Add", resid_old, changed_mol_counter, unchanged_mol_counter)
                        struct_calc.add_molecule(outf, resid_old, mol_pars)
                        added_mols_total += 1
                        resname_old = None
                else:
                    print(line, file=outf, end='')

main()
