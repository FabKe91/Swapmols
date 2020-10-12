'''
    Extracts information about system from given *.gro or *.pdb file.
'''
import logging
import numpy as np
import MDAnalysis as mda
from bilana import lipidmolecules
from .common import REGEX_PDB, REGEXP_GRO
LOGGER = logging.getLogger("Swapmols.systeminformation")

def count_mols_of(molname, struct_file_source):
    ''' Count how many molecules of type <molname> are in source structure (specified in -f) '''
    reached_once = False
    counter = 0
    counter_atoms =0
    resid_old = None
    with open(struct_file_source, "r") as topf:
        for line in topf:
            regmatch = REGEX_PDB.match(line)
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


def lipid_leaflet_assignment(grofilename, solvent='TIP3'):
    '''
        Reads grofile, searches all lipids and returns a dictionary
        assigning lipids to upper or lower leaflet
    '''
    leafleat_assign = {}
    coord_head = coord_base = None
    with open(grofilename, "r") as gfile:
        old_resid = 0
        sum_upper = 0
        sum_lower = 0
        for line in gfile:
            match = REGEXP_GRO.match(line)
            if match:
                resid = int(match.group(1).split()[0])
                if not old_resid:
                    old_resid = resid
                resname = match.group(2).split()[0].upper()
                atomname = match.group(3).split()[0]
                coords = [float(i) for i in match.group(5, 6, 7)]
                #logger.debug("Linepars: %s %s %s %s", resid, resname, atomname, coords)
                if resname == solvent:
                    continue
                last_tail_atm = lipidmolecules.scd_tail_atoms_of(resname)[0][-1]
                LOGGER.debug("Atmn/last_tail %s/%s", atomname, last_tail_atm)
                if old_resid != resid:
                    LOGGER.debug("Resid/resname/atmname %s/%s/%s", resid, resname, atomname)
                    LOGGER.debug("Coords head/base: %s/%s", coord_head, coord_base)
                    new_coords = coord_head - coord_base
                    cos = np.dot(new_coords, np.array([0.0,0.0,1.0]))/np.linalg.norm(new_coords)
                    if cos <= 0:
                        sum_upper += 1
                        leaflet = 0
                    else:
                        sum_lower += 1
                        leaflet = 1
                    leafleat_assign[old_resid] = leaflet
                    old_resid = resid
                    coord_head = coord_base = None
                if atomname in lipidmolecules.CENTRAL_ATOM_OF.values():
                    coord_head = np.array(coords)
                    LOGGER.debug("Added head %s", atomname)
                if atomname == last_tail_atm:
                    coord_base = np.array(coords)
                    LOGGER.debug("Added tail %s", atomname)
        if coord_base is not None and coord_head is not None:
            new_coords = coord_head - coord_base
            cos = np.dot(new_coords, np.array([0.0,0.0,1.0]))/np.linalg.norm(new_coords)
            if cos <= 0:
                sum_upper += 1
                leaflet = 0
            else:
                sum_lower += 1
                leaflet = 1
            leafleat_assign[old_resid] = leaflet
    return leafleat_assign


def radial_select(grofilename, mol_to_change, boxdivisor=3):
    '''
    Returns list of resids inside a radial area with radius box_x/boxdivisor around membrane center
    boxdivisor: Factor by which the box edge should be divided - default radius is one
                third of box.
    '''
    uni = mda.Universe(grofilename)
    resnames = set(uni.atoms.resnames)
    solvent_mols = []
    for solvent in ["TIP3", "SOL", "SPCE", "SPC", "TIP4", "NA", "K", "POT", "CL"]:
        try:
            resnames.remove(solvent)
            solvent_mols.append(solvent)
        except KeyError:
            continue
    LOGGER.debug("Solvent molecules %s", solvent_mols)
    if solvent_mols is None:
        solvent_mols = ''
    memb_sel = uni.select_atoms("resname {}".format(' '.join(resnames)))
    box = uni.dimensions
    membcenter = memb_sel.center_of_mass()
    LOGGER.debug("Box dimensions: %s", box)
    LOGGER.debug("Membrane center: %s", membcenter)

    #2 Get resids that are within r=thirdx around <membcenter>
    LOGGER.info("Getting list of resids within %s/%s of membrane center", box[0], boxdivisor)
    halfz = box[2]/2
    thirdx = box[0]/boxdivisor
    membcenter_str = "{} {} {}".format(*membcenter)
    LOGGER.debug("thirdx: %s, halfz %s, membcenter_str %s solvent_mols %s", thirdx, halfz, membcenter_str, solvent_mols)
    selection_str = "same residue as (cyzone {0} {1} -{1} (point {2} 10)) and resname {3}".format(thirdx, halfz, membcenter_str, mol_to_change)
    raft_atoms = uni.select_atoms(selection_str)
    raft_atoms.write("molecules_in_raft.gro")
    LOGGER.debug("Raft atom object: %s and number of resids %s", set(raft_atoms.resids), len(set(raft_atoms.resids)))
    return list(set(raft_atoms.resids))
