'''
    This module stores all things that use coordinate transformations to swap molecules
    Class CoordinateTransformation handles operations done with it

'''
import numpy as np
from src.logger import LOGGER
from src.common import REGEX_PDB, PDB_REC


class CoordinateTransformation():
    ''' Class for handling all coordinate transformation calculations
        functions:
            -add_molecule: Add all atoms from target structure to source structure
            using
                -gen_transmat: Calculate transformation matrix from target to source system
                -calc_coords: Calculates new coordinates for target atom names
            -get_target_top_information: Gets molecule name, atomname order, and dictionary
                                         atom name to coordinate from target structure

    '''
    def __init__(self,
                 targetstructure_fname,
                 refatmlist_src, refatmlist_targ,
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
        self.targetstructure_fname = targetstructure_fname

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
        returnvals = self.get_target_top_information(self.targetstructure_fname)
        self.atom_entry_order = returnvals[0]
        self.target_resname   = returnvals[1]
        self.atomname2relcoords_tar = returnvals[2]
        #
        self.number_of_atoms  = len(self.atom_entry_order)

    def add_molecule(self, fileobj, resid, atom_dict, other=''):
        ''' Main function to write all pdb lines for molecule <resid> in atom_entries variable
            refatom_dict: Dictionary for reference atoms forming coordinate system
                          linked to their coordinates
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
                regmatch = REGEX_PDB.match(line)
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



class Molecule():
    '''  '''
    def __init__(self):
        ''' '''
        self.finished = False
    def add_mol(self, fobj):
        ''' prints information to fobj '''
        self.finished = True
    def update_info(self, molinfo):
        ''' '''
        resid, resname = molinfo[1:3]
        atmname = molinfo[3] + molinfo[4]
        coords = molinfo[5:8]
