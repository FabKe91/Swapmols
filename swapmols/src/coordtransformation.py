'''
    This module stores all things that use coordinate transformations to swap molecules
    Class CoordinateTransformation handles operations done with it

'''
import logging
import numpy as np
from .common import REGEX_PDB, PDB_REC, GRO_STRING_FORMAT
from bilana import lipidmolecules

LOGGER = logging.getLogger("Swapmols.coordtransformation")

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
    '''
        Stores data of one residue that has to be changed
        and then prints changed mol to new structure file
    '''
    def __init__(self, swappair):
        '''
            swappair - list: [lip1, lip2] swap lip1 to lip2
            self.finished - determines wether stored information were used to add new molecule
        '''
        LOGGER.debug("Creating instance of molecule")
        self.finished = False
        self.swappair  = swappair
        self.missing_atoms_on_index = {}
        self.info_lines = []
        self.coord_dict = {}
        self.swapdict = self._get_swapdict(swappair)

    def add_mol(self,):
        ''' prints information to fobj '''
        LOGGER.debug("Adding molecule")
        LOGGER.debug("Swapdict: %s", self.swapdict)
        LOGGER.debug("missing atoms dict %s", self.missing_atoms_on_index)
        outputlines = []
        printafter = ''
        missingatm_printstring = ''
        for inf in self.info_lines:
            resid = int(inf[0])
            atmname = inf[2].strip()
            atomindex = inf[3]
            coords = np.array([float(i) for i in inf[4:7]])
            new_atmname = self._translate_atmname(atmname)
            if new_atmname is not None:
                outputlines.append(GRO_STRING_FORMAT.format(resid, self.swappair[1], new_atmname, atomindex,  *coords, ))

            c_serial = inf[2].strip()[1:]
            if c_serial in self.missing_atoms_on_index.keys():
                if missingatm_printstring:
                    print(missingatm_printstring)
                    raise ValueError("Would be overwriting above missing atom line")
                for missingatm in self.missing_atoms_on_index[c_serial]:
                    LOGGER.debug("missing atom: %s", missingatm)
                    if missingatm[0] == 'H':
                        if not printafter:
                            printafter = self._hydrogen_before(missingatm)
                        new_coords = self._add_hydrogen(atmname, coords)
                        missingatm_printstring += GRO_STRING_FORMAT\
                            .format(resid, self.swappair[1], missingatm, atomindex, *new_coords, )
                    else:
                        raise NotImplementedError("Can only add hydrogens.")
            if atmname == printafter:
                LOGGER.debug("Printafter %s,",  atmname)
                outputlines.append(missingatm_printstring)
                missingatm_printstring = printafter = ''
        self.finished = True
        return outputlines

    def update_info(self, molinfo):
        ''' Stores information given from stucture file
                self.coord_dict is basis for self._add_hydrogen()
                self.info_lines is looped through in self.add_mol()
        '''
        self.info_lines.append(molinfo)
        atmname = molinfo[2].strip()
        coords = np.array([float(i) for i in molinfo[4:7]])
        self.coord_dict[atmname] = coords

    def _get_swapdict(self, swappair):
        ''' '''
        def allequal(lst):
            return lst[1:] == lst[:-1]
        swapdict = {}
        is_sterol = {}
        sterol_in_pair = False
        for i, lip in enumerate(swappair):
            sterol = lipidmolecules.is_sterol(lip)
            is_sterol[i] = sterol
            if sterol:
                sterol_in_pair = True
        if sterol_in_pair:
            raise NotImplementedError("Not yet possible to switch phospholipid with sterols")
        else:
            tails = [swappair[0][:2], swappair[1][:2]]
            heads = [swappair[0][2:], swappair[1][2:]]
            LOGGER.debug("heads %s, tails %s", heads, tails)
            if not allequal(heads):
                headswap_dict = self._swapheads(swappair)
                swapdict.update(headswap_dict)
            if not allequal(tails):
                tailswap_dict = self._swapchains(swappair)
                swapdict.update(tailswap_dict)
        return swapdict


    def _swapchains(self, swappair):
        ''' Managing the swap of chains from swappair'''
        def convert_to_hdict(hydrogenlist):
            '''
                Converts hydrogen entries like H14S to {number:identifier}
                {14:S}
            '''
            chain_index = dict(R=2, S=2, T=2, X=3, Y=3, Z=3)
            converter = {}
            for hydr in hydrogenlist:
                index_carbon = hydr[1:-1]
                identifier = hydr[-1]
                number = str(chain_index[identifier])+index_carbon
                if number not in converter.keys():
                    converter[number] = []
                converter[number].append(hydr)
            return converter
        chaindict = {}
        tailcarbons_lip1 = lipidmolecules.tailcarbons_of(swappair[0])
        tailcarbons_lip2 = lipidmolecules.tailcarbons_of(swappair[1])
        tailhydr_lip1    = lipidmolecules.tailhydr_of(swappair[0])
        tailhydr_lip2    = lipidmolecules.tailhydr_of(swappair[1])
        if len(tailcarbons_lip1) == len(tailcarbons_lip2):# Check wether both lipids have two chains
            for chain, carbonlist_lip1 in enumerate(tailcarbons_lip1):
                hydr1_dict = convert_to_hdict(tailhydr_lip1[chain])
                hydr2_dict = convert_to_hdict(tailhydr_lip2[chain])
                if len(carbonlist_lip1) >= len(tailcarbons_lip2[chain]):# Loop over longer chain
                    for index, atom in enumerate(carbonlist_lip1):
                        skip_hydr = False
                        ci_chain =  atom[1:]
                        try: # Try, because atoms may be missing
                            cc_translation = {atom:tailcarbons_lip2[chain][index]}
                        except IndexError:
                            cc_translation = {atom:None}
                            skip_hydr = True
                        chaindict.update(cc_translation)
                        if not skip_hydr: # Dont look for missing H on carbons that are not added
                            hydrogens1_ci = hydr1_dict[ci_chain]
                            hydrogens2_ci = hydr2_dict[ci_chain]
                            if len(hydrogens1_ci) == len(hydrogens2_ci):
                                for h1, h2 in zip(hydrogens1_ci, hydrogens2_ci):
                                    hh_translation = {h1:h2}
                                    chaindict.update(hh_translation)
                            elif len(hydrogens1_ci) > len(hydrogens2_ci):
                                missing_h = set(hydrogens1_ci) ^ set(hydrogens2_ci)
                                self.missing_atoms_on_index[ci_chain] = missing_h
                                for hi, hydr in enumerate(hydrogens2_ci):
                                    hh_translation = {hydrogens1_ci[hi]:hydr}
                                    chaindict.update(hh_translation)
                            elif len(hydrogens1_ci) < len(hydrogens2_ci):
                                missing_h = set(hydrogens1_ci) ^ set(hydrogens2_ci)
                                self.missing_atoms_on_index[ci_chain] = missing_h
                                for hi, hydr in enumerate(hydrogens1_ci):
                                    hh_translation = {hydr:hydrogens2_ci[hi]}
                                    chaindict.update(hh_translation)
                        else:
                            hydrogens1_ci = hydr1_dict[ci_chain]
                            for hydr in hydrogens1_ci:
                                hh_translation = {hydr:None}
                                chaindict.update(hh_translation)
                else:
                    raise NotImplementedError("The target molecule's chain is longer"
                                              "than that of source molecule."
                                              "This is not yet implemented")

        return chaindict

    def _swapheads(self, swappair):
        raise NotImplementedError("Switching heads not yet implemented")

    def _add_hydrogen(self, atmn_c1, coords_c1):
        ''' Calculate new hydrogen coordinates based on tetrahedral geometry
                C0
             v1  |   \
                C02 -v2- C1 -v2- P +- |v1|*normal(C0-C1-C2)
                     /
                C2

        '''
        #tetrahedral_vectors = [
        #    np.array([ 1,  1,  1]),
        #    np.array([-1, -1,  1]),
        #    np.array([ 1, -1, -1]),
        #    np.array([-1,  1, -1]),
        #]

        # Get carbon before and after central carbon
        chain = atmn_c1[1]
        atmn_c0 = 'C{}{}'.format(chain, int(atmn_c1[2:])-1)
        atmn_c2 = 'C{}{}'.format(chain, int(atmn_c1[2:])+1)
        coords_c0 = self.coord_dict[atmn_c0]
        coords_c2 = self.coord_dict[atmn_c2]

        if self.swapdict[atmn_c2] is None:
            norm = np.linalg.norm(coords_c1-coords_c2)
            vec = (coords_c1-coords_c2)/norm * 0.11
            new_coords = coords_c1 - vec
            return new_coords

        # Calculate vectors v1, v2 and normal(C0-C1-C2)
        C02 = (coords_c0 + coords_c2)/2
        v1  = C02 - coords_c2
        v2  = coords_c1 - C02
        normal = np.cross(coords_c2-coords_c1, coords_c0-coords_c1)
        e_norm  = normal/np.linalg.norm(normal)
        LOGGER.debug("Parameters for calculation %s %s %s %s %s ", C02, v1, v2 , normal, e_norm)

        # Calculate positions of hydrogens
        for sign in [-1, 1]:
            h_coords = coords_c1+v2 + sign*np.linalg.norm(v1)/2 *e_norm
            too_close = False
            for key, coord in self.coord_dict.items():
                dist = np.linalg.norm(h_coords-coord)
                #print(key, coord, atmn_c1, dist)
                if dist <= 0.01:
                    too_close = True
            if not too_close:
                return h_coords
        raise ValueError("No new coordinate for hydrogen found.")

    def _translate_atmname(self, atmname):
        ''' translate name if in dict, else return atmname '''
        try:
            return self.swapdict[atmname]
        except KeyError:
            return atmname

    def _hydrogen_before(self, hydrogen):
        ''' Determines hydrogen atomname before <hydrogen>
            in charmm-ff structure files hydrogens are ordered like:
                HNX -> HNY -> HNZ
                HNR -> HNS -> HNT
        '''
        translator = dict(Z='Y', Y='X', T='S', S='R')
        identifier = hydrogen[-1]
        if identifier in ['X', 'R']:
            return '{}{}'.format('C', hydrogen[1:-1])
        else:
            return '{}{}'.format(hydrogen[:-1], translator[identifier])
