''' Reads user input flags and starts '''#### Parser arguments ####
import argparse

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
