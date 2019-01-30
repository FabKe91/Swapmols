'''
Add package description here

'''
import argparse
import numpy as np
import src

LOGGER = src.logger.LOGGER

PARSER = argparse.ArgumentParser()
PARSER.add_argument('-f', action='store', nARGS=1, metavar='input.pdb', required=True,
                    type=str, help='Input membrane structure file (pdb).')
PARSER.add_argument('-o', action='store', nARGS=1, metavar='output.pdb', required=True,
                    type=str, help='Output membrane structure file (pdb).')
PARSER.add_argument('-t', action='store', nARGS=1, metavar='structure_target.pdb', required=True,
                    type=str, help='Target molecule structure file (pdb).')
PARSER.add_argument('-i', action='store', nARGS=1, metavar='edit_specifications.inp', required=True,
                    type=str, help='File that contains a list of all atoms to be edited.')
PARSER.add_argument('-s', action='store', nARGS=1, metavar='<molname to modify>', required=True,
                    type=str, help='Name of molecule to change in source structure.\n'
                    'Syntax is described in function "read_input_atom_specs".')
PARSER.add_argument('-c', action='store', nARGS=1, metavar='<conc/%>', required=False,
                    type=int, help='Concentration of changed molecule per bilayer')
ARGS = PARSER.parse_args()   ### Arguments can be taken via args.<flagname>
INP_FNAME = ARGS.f[0]
OUT_FNAME = ARGS.o[0]
TARGSTRUCT_FNAME = ARGS.t[0]
EDITSPECS_FNAME = ARGS.i[0]
MOLNAME = ARGS.s[0]
ATMSEXPL = []
try:
    CONCENTRATION = ARGS.c[0]
except TypeError:
    CONCENTRATION = 100

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
    check_input()
    inp_specs = src.read_input_atom_specs(EDITSPECS_FNAME)
    name_mapping = {}
    changed_mol_counter = 0
    unchanged_mol_counter = 0
    added_mols_total = 0
    count_atoms = 0
    n_mols, n_atoms_mol = src.sysinfo.count_mols_of(MOLNAME, INP_FNAME)
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

    struct_calc = src.ctrans.CoordinateTransformation(
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

            regmatch = src.common.REGEX_PDB.match(line)
            if regmatch:
                grps = regmatch.groups()
                serial, atom_str, atom_int, resname, resid, x, y, z, other = grps
                atom = atom_str+atom_int    # letter and number of atomnames are handled differently

                # Start working on molecule or first molecule is already target:
                if resname_old == MOLNAME or (resname_old is None and resname == MOLNAME):
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
                if resname_old == MOLNAME and added_mols_total < n_mols:
                    if changed_mol_counter < mols_to_change:
                        LOGGER.info("Added %s mols and %s totals", n_mols, added_mols_total,)
                        LOGGER.info("Add %s - changed %s , unchanged %s",
                                    resid_old, changed_mol_counter, unchanged_mol_counter)
                        struct_calc.add_molecule(outf, resid_old, mol_pars)
                        added_mols_total += 1
                        resname_old = None
                else:
                    print(line, file=outf, end='')

def switch_lipidtails():
    ''' Transform lipid tails '''
    # Get all residues that are to be changed
    resids_to_change = src.sysinfo.radial_select(INP_FNAME)
    # Read whole gro file
    with open(INP_FNAME, "r") as fgro:
        grodata = fgro.readlines()
    # Morph residues in resid_to_change
    # otherwise print to OUT_FNAME
    with open(OUT_FNAME, "w") as outf:
        resid_old = None
        molecule = None
        for line in grodata:
            match =  src.common.REGEXP_GRO.match(line)
            if match:
                grps = match.groups()
                resid = grps[1]
                if resid_old != resid:# All information for resid_old gathered
                    if molecule is not None:# For first molecule
                        molecule.add_mol(outf)
                    molecule = src.coordtransformation.Molecule()
                if resid in resids_to_change:
                    molecule.update_info(grps)
                else:
                    if molecule is not None and not molecule.finished:
                        # If last residue not in resid to change
                        molecule.add_mol(outf)
                    # Print all lines until we get to new residue
                    print(line, file=outf)
            else:#Does not match if no atom entry
                # If last residue not in resid to change
                if molecule is not None and not molecule.finished:
                    molecule.add_mol(outf)
                print(line, file=outf)
