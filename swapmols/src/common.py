''' Generic functions and constants are stored here  '''
import re

#### GLOBALS ####
PDB_REC = 'ATOM'
#                    Record   Serial   atomname    resn    resid        x       y           z    Rest
PDB_STRING_FORMAT = '{: <4}  {: >5} {: >2}{: <2}{}{: >4}{}{: >4}{}   {: >8.3f}{: >8.3f}{: >8.3f}{}'#{: >6.2f}{: >6.2f}'
PDB_PATTERN = r"""
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
PDB_PATTERN = ''.join(PDB_PATTERN.split()) # Remove all whitespaces in above defined pattern
REGEX_PDB = re.compile(PDB_PATTERN)


#                     resid resnm atmnm atmnr   x        y         z         vx      vy        vz
GRO_STRING_FORMAT = '{: >5}{: <5}{: >5}{: >5}{: >8.3f}{: >8.3f}{: >8.3f}'#{: >8.4f}{: >8.4f}{: >8.4f} '
GRO_PATTERN = r"""
([\s,\d]{5})           (?# Resid       Grp 1)
([\s,\w]{5})           (?# Resname     Grp 2)
([\w,\s,\d,\']{5})           (?# Atom name   Grp 3)
([\s,\d]{5})           (?# Atom number Grp 4)
\s*
        (-?\d+\.\d+)      (?# X               Grp 5-7)
\s*
        (-?\d+\.\d+)      (?# Y)
\s*
        (-?\d+\.\d+)      (?# Z)
\s*
        (.*)            (?# The rest velocities Grp 8)
"""
GRO_PATTERN = ''.join(GRO_PATTERN.split())
REGEXP_GRO = re.compile(GRO_PATTERN)
GRO_BOX_PATTERN = r'\s*(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$'
REGEXP_BOX = re.compile(GRO_BOX_PATTERN)
