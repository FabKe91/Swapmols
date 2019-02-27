'''
Add package description here

'''
import argparse
import logging
import numpy as np
from .src import *

LOGGER = logging.getLogger("Swapmols.__main__")
SH, FH = mylogger.add_handlers(LOGGER)

PARSER = argparse.ArgumentParser()
PARSER.add_argument('-f', action='store', nargs=1, metavar='input.pdb', required=True,
                    type=str, help='Input membrane structure file (pdb).')
PARSER.add_argument('-o', action='store', nargs=1, metavar='output.pdb', required=True,
                    type=str, help='Output membrane structure file (pdb).')
PARSER.add_argument('-t', action='store', nargs=1, metavar='structure_target.pdb/lipidname', required=True,
                    type=str, help='Target molecule structure file (pdb) or name of lipid to switch to.')
PARSER.add_argument('-i', action='store', nargs=1, metavar='edit_specifications.inp', required=False,
                    type=str, help='File that contains a list of all atoms to be edited.'
                                    ' Only needed for mode molecule')
PARSER.add_argument('-s', action='store', nargs=1, metavar='<molname to modify>', required=True,
                    type=str, help='Name of molecule to change in source structure.\n'
                                   'Syntax is described in function "read_input_atom_specs".')
PARSER.add_argument('-c', action='store', nargs=1, metavar='<conc/%>', required=False,
                    type=int, help='Concentration of changed molecule per bilayer')
PARSER.add_argument('-m', action='store', nargs=1, metavar='lipid/molecule', required=True,
                    type=str, help='Sets mode: Either molecule changed using'
                                    'structure files or swap known lipids')
PARSER.add_argument('-v', action='count')
ARGS = PARSER.parse_args()   ### Arguments can be taken via args.<flagname>
MODE = ARGS.m[0]
INP_FNAME = ARGS.f[0]
OUT_FNAME = ARGS.o[0]
TARGSTRUCT_FNAME = ARGS.t[0]
MOLNAME_SRC = ARGS.s[0]
ATMSEXPL = []
if MODE == 'lipid' and ARGS.i is not None:
    LOGGER.warning("Warning: You gave edit specifications in -i with mode lipids."
                   "This flag is not needed in this mode.")
if MODE == 'molecule':
    EDITSPECS_FNAME = ARGS.i[0]
    try:
        CONCENTRATION = ARGS.c[0]
    except TypeError:
        CONCENTRATION = 100

mylogger.set_verbosity(ARGS.v, SH)


def check_input():
    ''' Checks if input is valid. Raise Error if not.'''
    if INP_FNAME == OUT_FNAME:
        raise ValueError("Input file must not be the same as output file:", INP_FNAME)
    elif TARGSTRUCT_FNAME == OUT_FNAME:
        raise ValueError("Input file must not be the same as output file:", TARGSTRUCT_FNAME)
    if INP_FNAME[-4:] != '.pdb' or TARGSTRUCT_FNAME[-4:] != '.pdb':
        print("WARNING are input structure files really in .pdb format?")


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

def switch_molecules():
    '''
        Reads input pdb file, gather data of molecules to morph (in mol_pars)
        and passes data to class transform_coordinates
    '''
    inp_specs = read_input_atom_specs(EDITSPECS_FNAME)
    name_mapping = {}
    changed_mol_counter = 0
    unchanged_mol_counter = 0
    added_mols_total = 0
    count_atoms = 0
    n_mols, n_atoms_mol = sysinfo.count_mols_of(MOLNAME_SRC, INP_FNAME)
    n_mols_per_bilayer = n_mols//2
    # In the case that only one molecule is to be changed, counter, counter_atoms
    if not n_mols_per_bilayer:
        n_mols_per_bilayer = 1
    mols_to_change = int(n_mols_per_bilayer*CONCENTRATION*0.01)
    LOGGER.info("Morphing %s molecules per bilayer of %s molecules to get a concentration of %s %%",
                mols_to_change, n_mols, CONCENTRATION)

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

    struct_calc = ctrans.CoordinateTransformation(
        TARGSTRUCT_FNAME,
        refatm_src, refatm_targ,
        skip_tar=skip_tar, mapping=name_mapping)

    LOGGER.info('Reference target: %s\nReference source:%sn'
                './%s\ntaking information from ./%s.\n\n',
                refatm_targ, refatm_src,
                INP_FNAME, TARGSTRUCT_FNAME)

    #### Start reading input and writing output bilayer structure ####
    with open(INP_FNAME, "r") as inpf, open(OUT_FNAME, "w") as outf:
        resid_old = 0     # To check if new residue sequence has begun
        resname_old = None
        first = True
        mol_pars = {}   # {<atmname>:XYZ}
        for line in inpf:

            regmatch = common.REGEX_PDB.match(line)
            if regmatch:
                grps = regmatch.groups()
                serial, atom_str, atom_int, resname, resid, x, y, z, other = grps
                atom = atom_str+atom_int    # letter and number of atomnames are handled differently

                # Start working on molecule or first molecule is already target:
                if resname_old == MOLNAME_SRC or (resname_old is None and resname == MOLNAME_SRC):
                    resid = int(resid)
                    xyz = np.array([float(val) for val in [x, y, z]])
                    mol_pars[atom] = xyz

                    # If enough molecules were morphed,
                    # but the residuals still have to be includeed:
                    if changed_mol_counter == mols_to_change and\
                            not (changed_mol_counter + unchanged_mol_counter == n_mols_per_bilayer):
                        if (count_atoms+1) / n_atoms_mol == 1:
                            LOGGER.info("Add source mol. Resid %s, Changed %s, Unchanged %s",
                                        resid_old, changed_mol_counter, unchanged_mol_counter)
                            unchanged_mol_counter += 1
                            count_atoms = 0
                        else:
                            print(line, file=outf, end='')
                            count_atoms += 1

                    # If all molecules were added
                    if changed_mol_counter + unchanged_mol_counter == n_mols_per_bilayer:
                        added_mols_total = added_mols_total\
                                           + changed_mol_counter\
                                           + unchanged_mol_counter
                        if unchanged_mol_counter != 0:# If not all residues are morphed
                            resid_old += 1
                        LOGGER.info("Reset. Nchanged %s and Nunchanged %s",
                                    changed_mol_counter, unchanged_mol_counter)
                        changed_mol_counter = 0
                        unchanged_mol_counter = 0
                        LOGGER.info("N total added: %s", added_mols_total)
                        #continue

                    # If new resid entry begun --> pass output.pdb-file object
                    # First var prevents script from immediately adding mol if
                    # struct file begins with mol
                    if resid_old != resid and not first:
                        if changed_mol_counter < mols_to_change and added_mols_total < n_mols:
                            LOGGER.info("Add %s - changed %s , unchanged %s",
                                        resid_old, changed_mol_counter, unchanged_mol_counter)
                            struct_calc.add_molecule(outf, resid_old, mol_pars)
                            changed_mol_counter += 1
                            mol_pars = {}

                    first = False
                    resname_old = resname

                else:
                    print(line, file=outf, end='')

                resid_old = resid
            else:
                 #Last calculation and no other residues come after
                if resname_old == MOLNAME_SRC and added_mols_total < n_mols:
                    if changed_mol_counter < mols_to_change:
                        LOGGER.info("Added %s mols and %s totals", n_mols, added_mols_total,)
                        LOGGER.info("Add %s - changed %s , unchanged %s",
                                    resid_old, changed_mol_counter, unchanged_mol_counter)
                        struct_calc.add_molecule(outf, resid_old, mol_pars)
                        added_mols_total += 1
                        resname_old = None
                else:
                    print(line, file=outf, end='')

def switch_lipids():
    '''
        Switches lipid types of input structure and write modified structure file
            --By now only switches of longer to shorter chains are implemented.
                * No switching of short to long lipid chain
                * No switching of heads

    '''
    new_lipid_lines = []
    swappair_resnames = [MOLNAME_SRC, TARGSTRUCT_FNAME]
    LOGGER.info("Switching %s to %s", MOLNAME_SRC, TARGSTRUCT_FNAME)
    # Get all residues that are to be changed
    resids_to_change = sysinfo.radial_select(INP_FNAME, MOLNAME_SRC)
    leaflet_assignment = sysinfo.lipid_leaflet_assignment("molecules_in_raft.gro")
    count = {0:0, 1:0}
    lipids_in_leaflet = {0:[], 1:[]}
    for key, val in leaflet_assignment.items():
        count[val] += 1
        lipids_in_leaflet[val].append(key)
    LOGGER.debug("Count: 0:%s  1:%s", count[0], count[1])
    if count[0] != count[1]:
        diff = count[0] - count[1]
        LOGGER.info("Bilayers not symmetric. Difference is %s", diff)
        LOGGER.info("Resids before: %s", resids_to_change)
        if diff < 0:
            del_resids = lipids_in_leaflet[1][diff:]
        elif diff > 0:
            del_resids = lipids_in_leaflet[0][diff*-1:]
        LOGGER.info("Deleting %s", del_resids)
        for i in del_resids:
            resids_to_change.remove(i)
        LOGGER.info("Resids after: %s", resids_to_change)
    # Read whole gro file
    with open(INP_FNAME, "r") as fgro:
        grodata = fgro.readlines()
    # Morph residues in resid_to_change
    # otherwise print to OUT_FNAME
    with open(OUT_FNAME, "w") as outf:
        resid_old = None
        resname_old = None
        molecule = None
        for line in grodata:
            match =  common.REGEXP_GRO.match(line)
            if match:
                grps = match.groups()
                resid = int(grps[0])
                resname = grps[1].strip()
                if resid_old != resid and resname == MOLNAME_SRC:# All information for resid_old gathered
                    resid_old = resid
                    if molecule is not None:# For first molecule
                        lipidlines = molecule.add_mol()
                        new_lipid_lines += lipidlines
                    molecule = ctrans.Molecule(swappair_resnames)
                if resname_old != resname:
                    print(resname_old)
                    if resname_old == MOLNAME_SRC:
                        LOGGER.info("Adding new molecules after %s and before %s", resname_old, resname)
                        for mol_entry in new_lipid_lines:
                            print(mol_entry, file=outf)
                    resname_old = resname
                    print(resname_old, "\n")
                if resid in resids_to_change:
                    molecule.update_info(grps)
                else:
                    #if molecule is not None and not molecule.finished:
                    #    # If last residue not in resid to change
                    #    molecule.add_mol(outf)
                    # Print all lines until we get to new residue
                    print(line, file=outf, end='')
            else:#Does not match if no atom entry
                # If last residue not in resid to change
                if molecule is not None:
                    LOGGER.debug("Adding last molecule")
                    lipidlines = molecule.add_mol()
                    if lipidlines:
                        new_lipid_lines.append(lipidlines)
                        for mol_entry in new_lipid_lines:
                            print(mol_entry, file=outf)

                print(line, file=outf, end='')

check_input()
if MODE == 'lipid':
    switch_lipids()
elif MODE == 'molecule':
    switch_molecules()
