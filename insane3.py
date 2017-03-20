#!/usr/bin/env python3

"""
INSANE: INSert membrANE
"""

import os
import sys
import math
import random
import collections
import argparse

version = "python3-beta"
previous = "20160625.15.TAW"

# ======================================================================
# 2.7 -> 3.5 conversion notes:
# Note: Code is now PEP8 compliant, ignoring line length limits.
# Note: Code contained multiple mirrored variables, these have
#       been converted from "var" to "_var" to avoid collisions.
# Note: Functions where documented when their function was clear.
# Note: File IO now uses the with-protocol.
# Note: main function was converted to class Insane3.
# Note: Option parsing has been implemented with argparse.
# Note: Unstable code has been marked as such.
# Note: Python 2's zip indexing has been fixed by doing, list(zip(a, b)) where needed
# Note: This build also addresses issues #1 and #5
# ======================================================================

# 20150920 - Insane can read lipid definitions from MARTINI repository:
#     ;@INSANE alhead=C P, allink=A A, altail=TCC CCCC, alname=DPSM, charge=0.0

# Set the random seed.
# The seed is set to an arbitary value set in the INSANE_SEED environment
# variable. If the environment variable is not set, then the system time is
# used to set the seed.
random.seed(os.environ.get('INSANE_SEED', None))

# Modify insane to take in arbitary lipid definition strings and use them as a template for lipids
# Also take in lipid name
# Edits: by Helgi I. Ingolfsson (all edits are marked with: # HII edit - lipid definition )

# PROTOLIPID (diacylglycerol), 18 beads
#
# 1-3-4-6-7--9-10-11-12-13-14
#  \| |/  |
#   2 5   8-15-16-17-18-19-20
#

lipidsx = {}
lipidsy = {}
lipidsz = {}
lipidsa = {}

#
# Diacyl glycerols
moltype = "lipid"
lipidsx[moltype] = (0, .5,  0,  0, .5,  0,  0, .5,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipidsy[moltype] = (0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipidsz[moltype] = (10,  9,  9,  8,  8,  7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)

lipidsa.update({
    #                   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
    # Phospholipids
    "DTPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    "DLPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DPPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DBPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "POPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    "DAPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 D1A D2A D3A D4A C5A  -  D1B D2B D3B D4B C5B  - "),
    "DIPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - "),
    "DGPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "DNPC": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B"),
    "DTPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    "DLPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DPPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DBPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "POPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPE": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    "POPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    "POPS": (moltype, " -   -   -  CNO  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B C2B C3B C4B  -   - "),
    "DOPS": (moltype, " -   -   -  CNO  -  PO4 GL1 GL2 C1A D2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    "DPSM": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A  -   -   -  C1B C2B C3B C4B  -   - "),
    "DBSM": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B  - "),
    "BNSM": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 T1A C2A C3A C4A  -   -  C1B C2B C3B C4B C5B C6B"),
    # PG for thylakoid membrane
    "OPPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B D2B C3B C4B  -   - "),
    # PG for thylakoid membrane of spinach (PPT with a trans-unsaturated bond at sn1 and a triple-unsaturated bond at sn2,
    # and PPG  with a transunsaturated bond at sn1 and a palmitoyl tail at sn2)
    "JPPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
    "JFPG": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  -   - "),
    # Monoacylglycerol
    "GMO":  (moltype, " -   -   -   -   -   -  GL1 GL2 C1A C2A D3A C4A C5A  -   -   -   -   -   -   - "),
    # Templates using the old lipid names and definitions
    "DHPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    "DMPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DSPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "POPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
    "DOPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "DUPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A D2A D3A C4A  -   -  C1B D2B D3B C4B  -   - "),
    "DEPC.o": (moltype, " -   -   -  NC3  -  PO4 GL1 GL2 C1A C2A C3A D4A C5A C6A C1B C2B C3B D4B C5B C6B"),
    "DHPE.o": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A  -   -   -   -  C1B C2B  -   -   -   - "),
    "DLPE.o": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DMPE.o": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A  -   -   -  C1B C2B C3B  -   -   - "),
    "DSPE.o": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "POPE.o": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
    "DOPE.o": (moltype, " -   -   -  NH3  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "PPCS.o": (moltype, " -   -   -  NC3  -  PO4 AM1 AM2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
    "DOPG.o": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "POPG.o": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
    "DOPS.o": (moltype, " -   -   -  CNO  -  PO4 GL1 GL2 C1A C2A D3A C4A C5A  -  C1B C2B D3B C4B C5B  - "),
    "POPS.o": (moltype, " -   -   -  CNO  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B C5B  - "),
    "CPG.o": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  C1B C2B D3B C4B  -   - "),
    "PPG.o": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A C2A C3A C4A  -   -  D1B C2B C3B C4B  -   - "),
    "PPT.o": (moltype, " -   -   -  GL0  -  PO4 GL1 GL2 C1A D2A D3A D4A  -   -  D1B C2B C3B C4B  -   - "),
    "DSMG.o": (moltype, " -   -   -  C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "DSDG.o": (moltype, "C61 C41 C11 C62 C42 C12 GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
    "DSSQ.o": (moltype, " -   -   S6 C6   C4 C1  GL1 GL2 C1A C2A C3A C4A C5A  -  C1B C2B C3B C4B C5B  - "),
})


# HII fix for PI templates and new templates PI(s) with diffrent tails, PO-PIP1(3) and POPIP2(4,5)
# Prototopology for phosphatidylinositol type lipids 5,6,7 are potentail phosphates (PIP1,PIP2 and PIP3)
# 1,2,3 - is the inositol and 4 is the phosphate that links to the tail part.
#  5
#   \
#  6-2-1-4-8--10-11-12-13-14-15
#    |/    |
#  7-3     9--16-17-18-19-20-21
moltype = "INOSITOLLIPIDS"
lipidsx[moltype] = (.5,  .5,   0,   0,   1, .5,  0,  0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipidsy[moltype] = (0,    0,   0,   0,   0,  0,  0,  0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipidsz[moltype] = (8,    9,   9,   7,  10, 10, 10,  6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipidsa.update({
    #                   1     2    3    4    5   6   7   8    9    10    11    12    13    14   15   16   17   18   19   20
    "OPI":  (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2   -    -    -    -    -    -   C1B  D2B  C3B  C4B   -    - "),
    "PPI":  (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2   -    -    -    -    -    -   C1B  C2B  C3B  C4B   -    - "),
    "DLPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  C3A   -    -    -   C1B  C2B  C3B   -    -    - "),
    "DPPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "DOPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  D2B  C3B  C4B   -    - "),
    "TPPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B   -    -    -    - "),
    "LPPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B   -    -    - "),
    "LOPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B   -    -    - "),
    "YPPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  C3B   -    -    - "),
    "POPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  D2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "PIPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  C1A  D2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "PAPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
    "PUPI": (moltype, " C1   C2   C3   PO4   -   -   -  GL1  GL2  D1A  D2A  D3A  D4A  D5A   -   C1B  C2B  C3B  C4B   -    - "),
    "POP1": (moltype, " C1   C2   C3   PO4  P1   -   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "POP2": (moltype, " C1   C2   C3   PO4  P1  P2   -  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "POP3": (moltype, " C1   C2   C3   PO4  P1  P2  P3  GL1  GL2  C1A  C2A  D3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # Templates using the old lipid names and definitions
    "PI.o":   (moltype, " C1   C2   C3    CP   -   -   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
    "PI34.o": (moltype, " C1   C2   C3    CP PO1 PO2   -  GL1  GL2  C1A  C2A  C3A  C4A   -    -   CU1  CU2  CU3  CU4  CU5   - "),
})


# Prototopology for IPC yeast lipid: MIP2C2OH, MIPC2OH, IPC2OH - Added by HII 2016.01.13
#
# 3-1-4-7-6
# |/     \|
# 2      5
#        |
#        9-8-11-12--14-15-16-17-18-19
#        \ |     |
#         10    13--20-21-22-23-24-25
moltype = "IPC"
lipidsx[moltype] = (1,   .5,   1,   1,  .5,   1,   1,  0,   0,  .5,    0,   0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipidsy[moltype] = (0,    0,   0,   0,   0,   0,   0,  0,   0,   0,    0,   0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipidsz[moltype] = (12,  13,  13,  11,   9,   9,  10,  8,   9,   8,    7,   6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipidsa.update({
    #                   1    2    3    4    5    6    7   8    9   10    11   12    13   14   15   16   17   18   19   20   21   22   23   24   25
    "PXI2": (moltype, "H31  H32  H33  H3P  H21  H22  H23 H11  H12  H13  PO4   AM1   AM2  O1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B  C5B  C6B"),
    "BXI2": (moltype, "H31  H32  H33  H3P  H21  H22  H23 H11  H12  H13  PO4   AM1   AM2  O1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B  C5B  C6B"),
    "PXI1": (moltype, " -    -    -    -   H21  H22  H23 H11  H12  H13  PO4   AM1   AM2  O1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B  C5B  C6B"),
    "PXI0": (moltype, " -    -    -    -    -    -    -  H11  H12  H13  PO4   AM1   AM2  O1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B  C5B  C6B"),
})


# Prototopology for longer and branched glycosil and ceramide based glycolipids
#
#     17-15-14-16
#         |/
#        13
#         |
# 12-10-9-7-6-4-3-1--18--20-21-22-23-24-25
#  |/   |/  |/  |/    |
#  11   8   5   2    19--26-27-28-29-30-31
moltype = "GLYCOLIPIDS"
lipidsx[moltype] = (0,  .5,   0,   0,  .5,  0,  0, .5,  0,    0,   .5,    0,    0,    0,   0,    0,    0,    0,   .5,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1)
lipidsy[moltype] = (0,   0,   0,   0,   0,  0,  0,  0,  0,    0,    0,    0,   .5,    1,   1,    1,    1,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0)
lipidsz[moltype] = (7,   8,   8,   9,  10, 10, 11, 12, 12,   13,   14,   14,   11,   10,  11,    9,   12,    6,    6,   5,   4,   3,   2,   1,   0,   5,   4,   3,   2,   1,   0)
lipidsa.update({
    #                   1     2    3    4    5   6   7   8   9    10    11    12    13    14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31
    "DPG1": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    "DXG1": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
    "PNG1": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    "XNG1": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
    "DPG3": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    "DXG3": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B  C5B  C6B"),
    "PNG3": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  D4B  C5B  C6B"),
    "XNG3": (moltype, "GM1  GM2  GM3  GM4  GM5 GM6  -   -   -    -     -     -    GM13  GM14 GM15  GM16  GM17   AM1   AM2  T1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  D4B  C5B  C6B"),
    "DPCE": (moltype, "  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    "DPGS": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  T1A  C2A  C3A   -    -    -   C1B  C2B  C3B  C4B   -    - "),
    "DPMG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "DPSG": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "DPGG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    # lipids for thylakoid membrane of cyanobacteria: oleoyl tail at sn1 and palmiotyl chain at sn2. SQDG no double bonds
    "OPMG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  C3B  C4B   -    - "),
    "OPSG": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  C3B  C4B   -    - "),
    "OPGG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  C3B  C4B   -    - "),
    # lipids for thylakoid membrane of spinach: for the *T both chains are triple unsaturated and the *G have a triple unsaturated chain at sn1 and a palmitoyl chain at sn2.
    "FPMG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  D3B  D4B   -    - "),
    "DFMG": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  D2A  D3A  D4A   -    -   C1B  D2B  D3B  D4B   -    - "),
    "FPSG": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  D3B  D4B   -    - "),
    "FPGG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  D2B  D3B  D4B   -    - "),
    "DFGG": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  D2A  D3A  D4A   -    -   C1B  D2B  D3B  D4B   -    - "),
    # Templates using the old lipid names and definitions
    "GM1.o":  (moltype, "GM1  GM2  GM3  GM4  GM5 GM6 GM7 GM8 GM9  GM10  GM11  GM12  GM13  GM14 GM15  GM16  GM17   AM1   AM2  C1A  C2A  C3A  C4A  C5A   -   C1B  C2B  C3B  C4B   -    - "),
    "DGDG.o": (moltype, "GB2  GB3  GB1  GA1  GA2 GA3   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "MGDG.o": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "SQDG.o": (moltype, " S1   C1   C2   C3    -   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "CER.o":  (moltype, "  -    -    -    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "GCER.o": (moltype, " C1   C2   C3    -    -   -   -   -   -     -     -     -     -     -    -     -     -   AM1   AM2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
    "DPPI.o": (moltype, " C1   C2   C3    -   CP   -   -   -   -     -     -     -     -     -    -     -     -   GL1   GL2  C1A  C2A  C3A  C4A   -    -   C1B  C2B  C3B  C4B   -    - "),
})


moltype = "QUINONES"
lipidsx[moltype] = (0,  .5,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipidsy[moltype] = (0,   0,   0,    0,   0,   0,   0,   0,   0,    0,    0,    0)
lipidsz[moltype] = (6,   7,   7,   5.5,  5,  4.5,  4,  3.5, 2.5,   2,  1.5,    1)
lipidsa.update({
    #                   1     2    3    4    5    6    7    8    9    10    11    12
    "PLQ": (moltype, " PLQ3 PLQ2 PLQ1 PLQ4 PLQ5 PLQ6 PLQ7 PLQ8 PLQ9 PLQ10 PLQ11 PLQ12"),
})


# Prototopology for cardiolipins
#
#       4-11-12-13-14-15-16
#       |
#   2---3--5--6--7--8--9-10
#  /
# 1
#  \
#   17-18-20-21-22-23-24-25
#       |
#      19-26-27-28-29-30-31
#
moltype = "CARDIOLIPINS"
lipidsx[moltype] = (0.5,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1)
lipidsy[moltype] = (1,     0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,   0,  0,  1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1)
lipidsz[moltype] = (8,     7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0,   7,  6,  6,  5,  4,  3,  2,  1,  0,  5,  4,  3,  2,  1,  0)
lipidsa.update({
    #                   1    2   3   4   5   6   7   8   9  10  11  12  13  14  15  16   17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
    "CDL0": (moltype, "GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), # Warning not the same names is in .itp
    "CDL1": (moltype, "GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), # Warning not the same names is in .itp
    "CDL2": (moltype, "GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"), # Warning not the same names is in .itp
    "CL4P": (moltype, "GL5 PO41 GL1 GL2 C1A C2A C3A C4A C5A   - C1B C2B C3B C4B C5B   - PO42 GL3 GL4 C1C C2C C3C C4C C5C   - C1D C2D C3D C4D C5D   -"),
    "CL4M": (moltype, "GL5 PO41 GL1 GL2 C1A C2A C3A   -   -   - C1B C2B C3B   -   -   - PO42 GL3 GL4 C1C C2C C3C   -   -   - C1D C2D C3D   -   -   -"),
    # Templates using the old lipid names and definitions
    "CL4.o":  (moltype, "GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"),
    "CL4O.o": (moltype, "GL5 PO41 GL1 GL2 C1A C2A D3A C4A C5A   - C1B C2B D3B C4B C5B   - PO42 GL3 GL4 C1C C2C D3C C4C C5C   - C1D C2D D3D C4D C5D   -"),
})


# Prototopology for mycolic acid(s)
#
#  1--2--3--4--5--6--7--8
#                       |
# 16-15-14-13-12-11-10--9
# |
# 17-18-19-20-21-22-23-24
#                     /
# 32-31-30-29-28-27-25-26
#

moltype = "MYCOLIC ACIDS"
lipidsx[moltype] = (0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    1,    1,    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1)
lipidsy[moltype] = (0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,    1,    1,    1,    1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,  0,   0,   0)
lipidsz[moltype] = (7,   6,   5,   4,   3,   2,   1,   0,   0,   1,   2,   3,   4,   5,   6,    7,    7,    6,    5,   4,   3,   2,   1,   0,   1,   0,   2,   3,   4,  5,   6,   7)
lipidsa.update({
    #                     1    2    3    4    5    6    7    8    9   10   11   12   13   14   15    16    17    18    19   20   21   22   23   24   25   26   27   28   29   30   31   32
    "AMA":   (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    "AMA.w": (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    "KMA":   (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
    "MMA":   (moltype, "  -    -    -  C1A  C2A  C3A  C4A  C5A  M1A  C1B  C2B  C3B  C4B    -    -     -     -     -   M1B  C1C  C2C  C3C    -    -  COH  OOH  C1D  C2D  C3D  C4D  C5D  C6D"),
})


# Sterols
moltype = "sterol"
lipidsx[moltype] = (0,  0,  0,  0,  0, 0,   0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0)
lipidsy[moltype] = (0,  0,  0,  0,  0, 0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
lipidsz[moltype] = (0,  0,  0,  0,  0, 0, 5.3,4.5,3.9,3.3, 3 ,2.6,1.4,  0,  0,  0,  0,  0)
lipidsa.update({
    "CHOL": (moltype, " -   -   -   -   -   -  ROH  R1  R2  R3  R4  R5  C1  C2  -   -   -   - "),
    "ERGO": (moltype, " -   -   -   -   -   -  ROH  R1  R2  R3  R4  R5  C1  C2  -   -   -   - "),
})


# Hopanoids
moltype = "Hopanoids"
lipidsx[moltype] = (0,  0,  0,  0, 0.5, -0.5,   0,   0, 0.5, 0.5,   0,   0,   0,   0,  0,  0,  0,  0)
lipidsy[moltype] = (0,  0,  0,  0,   0,    0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0)
lipidsz[moltype] = (0,  0,  0,  0, 0.5,  1.4, 2.6,   3, 3.3, 3.9, 4.5, 5.0, 5.5, 6.0,  0,  0,  0,  0)
lipidsa.update({
    "HOPR": (moltype, " -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   -    -    -    -   -   -   - "),
    "HHOP": (moltype, " -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   -    -    -   -   -   - "),
    "HDPT": (moltype, " -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   -    -    -   -   -   - "),
    "HBHT": (moltype, " -   -   -   R1   R2   R3   R4   R5   R6   R7   R8   C1   C2   C3   -   -   -   - "),
})


# Lists for automatic charge determination
charges = {
    "ARG": 1, "LYS": 1, "ASP": -1, "GLU": -1, "DOPG": -1,
    "POPG": -1, "DOPS": -1, "POPS": -1, "DSSQ": -1
}

a,  b = math.sqrt(2) / 20, math.sqrt(2) / 60
ct, st = math.cos(math.pi * 109.47 / 180), math.sin(math.pi * 109.47 / 180)  # Tetrahedral

# Get a set of coordinates for a solvent particle with a given name
# Dictionary of solvents; First only those with multiple atoms
solventParticles = {
    "PW":       (("W", (-0.07, 0, 0)),                          # Polarizable water
                 ("WP", (0.07, 0, 0)),
                 ("WM", (0.07, 0, 0))),
    "BMW":      (("C", (0, 0, 0)),
                 ("Q1", (0.12, 0, 0)),
                 ("Q2", (-0.06, math.cos(math.pi / 6) * 0.12, 0))),  # BMW water
    "SPC":      (("OW", (0, 0, 0)),                             # SPC
                 ("HW1", (0.01, 0, 0)),
                 ("HW2", (0.01 * ct, 0.01 * st, 0))),
    "SPCM":     (("OW", (0, 0, 0)),                             # Multiscale/Martini SPC
                 ("HW1", (0.01, 0, 0)),
                 ("HW2", (0.01 * ct, 0.01 * st, 0)),
                 ("vW", (0, 0, 0))),
    "FG4W": (("OW1", (a, a, a)),  # Bundled water
             ("HW11", (a, a - b, a - b)),
             ("HW12", (a, a + b, a + b)),
             ("OW2", (a, -a, -a)),
             ("HW21", (a - b, -a, -a + b)),
             ("HW22", (a + b, -a, -a - b)),
             ("OW3", (-a, -a, a)),
             ("HW31", (-a, -a + b, a - b)),
             ("HW32", (-a, -a - b, a + b)),
             ("OW4", (-a, a, -a)),
             ("HW41", (-a + b, a, -a + b)),
             ("HW42", (-a - b, a, -a - b))),
    "FG4W-MS": (("OW1", (a, a, a)),  # Bundled water, multiscaled
                ("HW11", (a, a - b, a - b)),
                ("HW12", (a, a + b, a + b)),
                ("OW2", (a, -a, -a)),
                ("HW21", (a - b, -a, -a + b)),
                ("HW22", (a + b, -a, -a - b)),
                ("OW3", (-a, -a, a)),
                ("HW31", (-a, -a + b, a - b)),
                ("HW32", (-a, -a - b, a + b)),
                ("OW4", (-a, a, -a)),
                ("HW41", (-a + b, a, -a + b)),
                ("HW42", (-a - b, a, -a - b)),
                ("VZ", (0, 0, 0))),
    "GLUC": (("B1", (-0.11, 0, 0)),
             ("B2", (0.05, 0.16, 0)),
             ("B3", (0.05, -0.16, 0))),
    "FRUC": (("B1", (-0.11, 0, 0)),
             ("B2", (0.05, 0.16, 0)),
             ("B3", (0.05, -0.16, 0))),
    "SUCR": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "MALT": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "CELL": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "KOJI": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "SOPH": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "NIGE": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "LAMI": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    "TREH": (("B1", (-0.25, 0.25, 0)),
             ("B2", (-0.25, 0, 0)),
             ("B3", (-0.25, -0.25, 0)),
             ("B4", (0.25, 0, 0)),
             ("B5", (0.25, 0.25, 0)),
             ("B6", (0.25, -0.25, 0))),
    # Loose aminoacids
    "GLY": (("BB", (0, 0, 0)),),
    "ALA": (("BB", (0, 0, 0)),),
    "ASN": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "ASP": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "GLU": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "GLN": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "LEU": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "ILE": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "VAL": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "SER": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "THR": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "CYS": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "MET": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "LYS": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "PRO": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "HYP": (("BB", (0.25, 0, 0)),
            ("SC1", (-0.25, 0, 0))),
    "ARG": (("BB", (0.25, 0, 0)),
            ("SC1", (0, 0, 0)),
            ("SC2", (-0.25, 0.125, 0))),
    "PHE": (("BB", (0.25, 0, 0)),
            ("SC1", (0, 0, 0)),
            ("SC2", (-0.25, -0.125, 0)),
            ("SC3", (-0.25, 0.125, 0))),
    "TYR": (("BB", (0.25, 0, 0)),
            ("SC1", (0, 0, 0)),
            ("SC2", (-0.25, -0.125, 0)),
            ("SC3", (-0.25, 0.125, 0))),
    "TRP": (("BB", (0.25, 0.125, 0)),
            ("SC1", (0.25, 0, 0)),
            ("SC2", (0, -0.125, 0)),
            ("SC3", (0, 0.125, 0)),
            ("SC4", (-0.25, 0, 0))),
    }

# Update the solvents dictionary with single atom ones
for s in ("W", "NA", "CL", "Mg", "K", "BUT"):
    solventParticles[s] = ((s, (0, 0, 0)),)

# Apolar amino acids nd stuff for orienting proteins in membrane
apolar = "ALA CYS PHE ILE LEU MET VAL TRP PLM CLR".split()

# PRIVATE PARTS FROM THIS POINT ON ##

S = str
F = float
I = int
R = random.random


def vector(v):
    """
    Convert string to floating point vector.
    :param v: string vector e.g. "1.23" or "2.5,1.6,-8.74"
    :return: Array<float> or float
    """
    if type(v) == str and "," in v:
        return [float(i) for i in v.split(",")]
    return float(v)


def vvadd(_a, _b):
    """
    Add two vectors
    :param _a: vector A
    :param _b: vector B
    :return: Array<float>
    """
    if type(_b) in (int, float):
        return [i+_b for i in _a]
    return [i+j for i, j in zip(_a, _b)]


def vvsub(_a, _b):
    """
    Subtract two vectors
    :param _a: vector A
    :param _b: vector B
    :return: Array<float>
    """
    if type(_b) in (int, float):
        return [i-_b for i in _a]
    return [i-j for i, j in zip(_a, _b)]


def mean(_a):
    """
    Calculate the mean.
    :param _a: vector of size n
    :return: float
    """
    return sum(_a) / len(_a)


def is_pdb_atom(l):  # Name needs to be refactored
    """
    Check if line represents an atom.
    :param l: line
    :return: boolean
    """
    return l.startswith("ATOM") or l.startswith("HETATM")


def pdb_atom(a):  # Name needs to be refactored
    """
    Convert a PDB atom line into useful objects
    :param a: atom line
    :return: Array<String|float|int>
    """
    # 01234567890123456789012345678901234567890123456789012345678901234567890123456789
    # ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    #  ===>   atom name,   res name,     res id, chain,       x,            y,             z
    return (S(a[12:16]), S(a[17:20]), I(a[22:26]), a[21], F(a[30:38]) / 10,
            F(a[38:46]) / 10, F(a[46:54]) / 10)


# Reformatting of lines in structure file
pdbBoxLine = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"
pdbline = "ATOM  %5i  %-3s %4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s  \n"

d2r = 3.14159265358979323846264338327950288 / 180


def norm2(_a):
    """
    Calculate the sum of squares for a vector.
    :param _a: vector
    :return: float
    """
    return sum([i*i for i in _a])


def norm(_a):
    """
    Calculate the square root of the square sum for a vector.
    :param _a: vector
    :return: float
    """
    return math.sqrt(norm2(_a))


def cos_angle(_a, _b):
    """
    Calculate cosine-angle between two vectors.
    :param _a: vector
    :param _b: vector
    :return: float
    """
    p = sum([i*j for i, j in zip(_a, _b)])
    q = math.sqrt(sum([i*i for i in _a]) * sum([j*j for j in _b]))
    return min(max(-1, p / q), 1)


def pdb_box_string(box):  # Name needs to be refactored
    """
    Calculate the PDB Crystal dimensions and format string.
    :param box: vector of n=3
    :return: String
    """
    # Box vectors
    u, v, w = box

    # Box vector lengths
    nu, nv, nw = [math.sqrt(norm2(i)) for i in (u, v, w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(cos_angle(v, w)) / d2r
    beta = nu*nw == 0 and 90 or math.acos(cos_angle(u, w)) / d2r
    gamma = nu*nv == 0 and 90 or math.acos(cos_angle(u, v)) / d2r

    return pdbBoxLine % (10 * norm(u), 10 * norm(v), 10 * norm(w), alpha, beta, gamma)


def pdb_box_read(_a):  # Name needs to be refactored
    """
    Convert a PDB CRYST1 entry to a lattice definition.
    Convert from Angstrom to nanometer
    :param _a: vector
    :return: Array<float|int>
    """
    fa, fb, fc, aa, ab, ac = [float(i) for i in _a.split()[1:7]]
    ca, cb, cg, sg = math.cos(d2r * aa), math.cos(d2r * ab), math.cos(d2r * ac), math.sin(d2r * ac)
    wx, wy = 0.1 * fc * cb, 0.1 * fc * (ca - cb * cg) / sg
    wz = math.sqrt(0.01 * fc * fc - wx * wx - wy * wy)
    return [0.1 * fa, 0, 0, 0.1 * fb * cg, 0.1 * fb * sg, 0, wx, wy, wz]


def gro_atom(_a):  # Name needs to be refactored
    """
    Convert a GRO atom line into useful objects
    :param _a: line
    :return: Array<String|float|int>
    """
    # 012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    #  ===>   atom name,   res name,     res id, chain,       x,          y,          z
    return S(_a[10:15]), S(_a[5:10]), I(_a[:5]), " ", F(_a[20:28]), F(_a[28:36]), F(_a[36:44])


def gro_box_read(_a):  # Name needs to be refactored
    """
    Needs description.
    :param _a: string
    :return: Array<float>
    """
    _b = [F(i) for i in _a.split()] + 6 * [0]  # Padding for rectangular boxes
    return _b[0], _b[3], _b[4], _b[5], _b[1], _b[6], _b[7], _b[8], _b[2]


def read_box(_a):  # Name needs to be refactored
    """
    Get lattice definition.
    :param _a: string
    :return: Array<float>
    """
    x = [float(i) for i in _a.split(",")] + 6 * [0]
    if len(x) == 12:  # PDB format
        return pdb_box_read("CRYST1 " + " ".join([str(i) for i in x]))
    else:  # GRO format
        return x[0], x[3], x[4], x[5], x[1], x[6], x[7], x[8], x[2]


class Structure:

    def __init__(self, filename=None):
        self.title = ""
        self.atoms = []
        self.coord = []
        self.rest = []
        self.box = []
        self._center = None

        # Note: This should probably be its own method.
        if filename:
            lines = []
            with open(filename) as _fh:
                lines = _fh.readlines()
            # Try extracting PDB atom/hetatm definitions
            self.rest = []
            self.atoms = [pdb_atom(i) for i in lines if is_pdb_atom(i) or self.rest.append(i)]
            if self.atoms:
                # This must be a PDB file
                self.title = "THIS IS INSANE!\n"
                for i in self.rest:
                    if i.startswith("TITLE"):
                        self.title = i
                self.box = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                for i in self.rest:
                    if i.startswith("CRYST1"):
                        self.box = pdb_box_read(i)
            else:
                # This should be a GRO file
                self.atoms = [gro_atom(i) for i in lines[2:-1]]
                self.rest = [lines[0], lines[1], lines[-1]]
                self.box = gro_box_read(lines[-1])
                self.title = lines[0]
            self.coord = [i[4:7] for i in self.atoms]
            self.center()

    def __nonzero__(self):
        return bool(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __iadd__(self, s):
        for i in range(len(self)):
            self.coord[i] = vvadd(self.coord[i], s)
        return self

    def center(self, other=None):
        if not self._center:
            self._center = [sum(i) / len(i) for i in zip(*self.coord)]
        if other:
            s = vvsub(other, self._center)
            for i in range(len(self)):
                self.coord[i] = vvadd(self.coord[i], s)
            self._center = other
            return s  # return the shift
        return self._center

    def diam(self):
        if self._center != (0, 0, 0):
            self.center((0, 0, 0))
        return 2 * math.sqrt(max([i * i + j * j + k * k for i, j, k in self.coord]))

    def diamxy(self):
        if self._center != (0, 0, 0):
            self.center((0, 0, 0))
        return 2 * math.sqrt(max([i * i + j * j for i, j, k in self.coord]))

    def fun(self, fn):
        return [fn(i) for i in zip(*self.coord)]


headbeads = {
    # Define supported lipid head beads. One letter name mapped to atom name
    "C":  "NC3",  # NC3 = Choline
    "E":  "NH3",  # NH3 = Ethanolamine
    "G":  "GL0",  # GL0 = Glycerol
    "S":  "CNO",  # CNO = Serine
    "P":  "PO4",  # PO4 = Phosphate
    }
linkbeads = {
    # Define supported lipid link beads. One letter name mapped to atom name
    "G":  "GL",  # Glycerol
    "A":  "AM",  # Amide (Ceramide/Sphingomyelin)
}


class Lipid:
    """Lipid structure"""

    def __init__(self, **kwargs):
        self.name = kwargs.get("name")
        self.head = kwargs.get("head")
        self.link = kwargs.get("link")
        self.tail = kwargs.get("tail")
        self.beads = kwargs.get("beads")
        if type(self.beads) == str:
            self.beads = self.beads.split()
        self.charge = kwargs.get("charge")
        self.template = kwargs.get("template")
        self.area = kwargs.get("area")
        self.diam = kwargs.get("diam", math.sqrt(kwargs.get("area", 0)))
        self.coords = None
        if kwargs.get("string"):
            self.parse(kwargs["string"])

    def parse(self, string):
        """
        Parse lipid definition from string:

            alhead=C P, allink=A A, altail=TCC CCCC, alname=DPSM, charge=0.0
        """
        fields = [i.split("=") for i in string.split(',')]
        for what, val in fields:
            what = what.strip()
            val = val.split()
            if what.endswith("head"):
                self.head = val
            elif what.endswith("link"):
                self.link = val
            elif what.endswith("tail"):
                self.tail = val
            elif what == "charge":
                self.charge = float(val[0])
            elif what.endswith("name") and not self.name:
                self.name = val[0]
        if self.charge is None:
            # Infer charge from head groups
            # Note: This code will crash, headgroup_charges does not exist
            self.charge = sum([headgroup_charges[bead] for bead in self.head])

    def build(self, **kwargs):
        """Build/return a list of [(bead,x,y,z),...]"""

        if not self.coords:
            if self.beads and self.template:
                stuff = zip(self.beads, self.template)
                self.coords = [[i, x, y, z] for i, (x, y, z) in stuff if i != "-"]
            else:
                # Set beads/structure from head/link/tail
                # Set bead names
                if self.beads:
                    beads = list(self.beads)
                else:
                    beads = [headbeads[i] for i in self.head]
                    beads.extend([linkbeads[n] + str(i + 1) for i, n in enumerate(self.link)])
                    for i, t in enumerate(self.tail):
                        beads.extend([n + chr(65 + i) + str(j + 1) for j, n in enumerate(t)])

                taillength = max([0] + [len(i) for i in self.tail])
                length = len(self.head) + taillength

                # Add the pseudocoordinates for the head
                rl = range(len(self.head))
                struc = [(0, 0, length - i) for i in rl]

                # Add the linkers
                rl = range(len(self.link))
                struc.extend([(i % 2, i / 2, taillength) for i in rl])

                # Add the tails
                for j, tail in enumerate(self.tail):
                    rl = range(len(tail))
                    struc.extend([(j % 2, j / 2, taillength - 1 - i) for i in rl])

                mx, my, mz = [(max(i) + min(i)) / 2 for i in zip(*struc)]
                self.coords = [[i, 0.25 * (x - mx), 0.25 * (y - my), z] for i, (x, y, z) in list(zip(beads, struc))]

        # Scale the x/y based on the lipid's APL - diameter is less than sqrt(APL)
        diam = kwargs.get("diam", self.diam)
        radius = diam * 0.45
        minmax = [(min(i), max(i)) for i in list(zip(*self.coords))[1:]]
        mx, my, mz = [sum(i) / 2. for i in minmax]
        scale = radius / math.sqrt((minmax[0][0] - mx) ** 2 + (minmax[1][0] - my) ** 2)

        for i in self.coords:
            i[1] = scale * (i[1] - mx)
            i[2] = scale * (i[2] - my)
            i[3] -= minmax[2][0]

        return self.coords

    def h(self, head):
        self.head = head.replace(".", " ").split()

    def l(self, link):
        self.link = link.replace(".", " ").split()

    def t(self, tail):
        self.tail = tail.replace(".", " ").split()

    def c(self, charge):
        self.charge = float(charge)


class LipidList(collections.MutableMapping):  # Name needs to be refactored
    """Container class for lipid definitions"""

    def __init__(self):
        self.store = dict()
        self.last = None

    def __getitem__(self,key):
        if key == -1:
            return self.last
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def add(self, name=None, string=None):
        lip = Lipid(name=name, string=string)
        self.store[lip.name] = lip
        self.last = lip.name


def meand(v):
    """
    Mean of deviations from initial value.
    :param v: initial value
    :return: float
    """
    return sum([i-v[0] for i in v]) / len(v)


def ssd(u, v):
    """
    Sum of squares / cross-products of deviations.
    :param u: vector U
    :param v: vector V
    :return: float
    """
    return sum([(i - u[0]) * (j - v[0]) for i, j in zip(u, v)]) / (len(u) - 1)


def parse_mol(x):
    """
    Parse a string for a lipid as given on the command line (LIPID[:NUMBER]).
    :param x: string of form LIPID[:NUMBER]
    :return: String and float|int
    """
    l = x.split(":")
    return l[0], len(l) == 1 and 1 or float(l[1])


# === MIJN EIGEN ROUTINE ===

# Quite short piece of code for diagonalizing symmetric 3x3 matrices :)

def solve_p3(_a, _b, _c):
    """
    Analytic solution for third order polynomial.
    :param _a: vector A
    :param _b: vector B
    :param _c: vector C
    :return: Array<float>
    """
    q, _r, a3 = (3 * _b - _a ** 2) / 9.0, (-27 * _c + _a * (9 * _b - 2 * _a ** 2)) / 54.0, _a / 3.0
    if q ** 3 + _r ** 2:
        t, r13 = math.acos(_r / math.sqrt(-q ** 3)) / 3, 2 * math.sqrt(-q)
        u, v, w = math.cos(t), math.sin(t + math.pi / 6), math.cos(t + math.pi / 3)
        return r13 * u - a3, -r13 * v - a3, -r13 * w - a3
    else:
        # Note: This code will crash! No such method sqrt3
        r13 = math.sqrt3(_r)
        return 2 * r13 - a3, -r13 - a3, -r13 - a3


def normalize(_a):
    """
    Normalization of 3-vector
    :param _a: vector of size n=3
    :return: Array<float>
    """
    f = 1.0 / math.sqrt(_a[0] * _a[0] + _a[1] * _a[1] + _a[2] * _a[2])
    return f * _a[0], f * _a[1], f * _a[2]


def mijn_eigen_sym_3x3(_a, d, f, _b, c, e):
    """
    Eigenvectors for a symmetric 3x3 matrix:
    For symmetric matrix A the eigenvector v with root r satisfies
    v.Aw = Av.w = rv.w = v.rw
    v.(A-rI)w = v.Aw - v.rw = 0 for all w
    This means that for any two vectors p,q the eigenvector v follows from:
    (A-rI)p x (A-rI)q
    The input is var(x),var(y),var(z),cov(x,y),cov(x,z),cov(y,z)
    The routine has been checked and yields proper eigenvalues/-vectors
    :param _a: ?
    :param d: ?
    :param f: ?
    :param _b: ?
    :param c: ?
    :param e: ?
    :return: ?
    """
    _a, d, f, _b, c, e = 1, d / _a, f / _a, _b / _a, c / _a, e / _a
    b2, c2, e2, df = _b * _b, c * c, e * e, d * f
    roots = list(
        solve_p3(-_a - d - f, df - b2 - c2 - e2 + _a * (f + d), _a * e2 + d * c2 + f * b2 - _a * df - 2 * _b * c * e))
    roots.sort(reverse=True)
    ux, uy, uz = _b * e - c * d, _b * c - _a * e, _a * d - _b * _b
    u = (ux + roots[0] * c, uy + roots[0] * e, uz + roots[0] * (roots[0] - _a - d))
    v = (ux + roots[1] * c, uy + roots[1] * e, uz + roots[1] * (roots[1] - _a - d))
    w = u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]  # Cross product
    return normalize(u), normalize(v), normalize(w), roots


# This class replaces the crowded main function from previous versions
class Insane3:

    def __init__(self, argv=None):
        # Variables
        self.options = {}
        # Lipid variables
        self.tm = []
        self.lipL, self.lipU, self.solv = [], [], []
        self.usrlip = LipidList()
        self.mollist, self.usrnames, self.usrheads = [], [], []
        self.usrlinks, self.usrtails, self.usrcharg = [], [], []
        self.parse()
        # PBC variables
        self.pbc_set_x, self.pbc_set_y, self.pbc_set_z = 0, 0, 0
        self.pbcx, self.pbcy, self.pbcz = 0, 0, 0
        self.lo_lipd = math.sqrt(self.options["a"])
        self.up_lipd = math.sqrt(self.options["au"]) if self.options["au"] else self.lo_lipd
        # Protein variables
        self.protein = Structure()
        self.prot = []
        self.xshifts = [0]
        self.prot_up, self.prot_lo = [], []
        self.resi = 0
        self.atid, self.molecules, self.box = 0, [], []
        self.grobox, self.rx, self.ry, self.rz = [], 0, 0, 0
        # Membrane variables
        self.membrane = Structure()
        self.numU, self.numL = 0, 0
        # Solvent variables
        self.solvent = None
        self.sol = []
        # Export variables
        self.title = None
        # Protocol
        self.tm = [Structure(t) for t in self.tm]
        self.absoluteNumbers = not self.options["d"]
        self.liplist = self.process_lipids()
        self.process_box()
        self.process_protein()
        self.process_membrane()
        self.process_solvents()
        self.write_to_gromacs_or_pdb()
        self.topology = self.__get_topology()
        if self.options["p"]:
            self.write_topology()
        else:
            print("\n".join("%-10s %7d"%i for i in self.topology), file=sys.stderr)

    def parse(self):
        """
        Add options to the command line.
        :return: None
        """
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("-v", "--version", action="version", version=version)
        # Add IO arguments
        parser.add_argument("-f", required=False, nargs="+", type=str, default=None, action="append", help="Input GRO or PDB file 1: Protein")  # Action
        parser.add_argument("-o", required=True, type=str, default=None, help="Output GRO file: Membrane with Protein")
        parser.add_argument("-p", required=False, type=str, default=None, help="Optional rudimentary topology file")
        # Periodic boundary conditions
        """
            Periodic boundary conditions
            If -d is given, set up PBC according to -pbc such that no periodic
            images are closer than the value given.  This will make the numbers
            provided for lipids be interpreted as relative numbers. If -d is
            omitted, those numbers are interpreted as absolute numbers, and the
            PBC are set to fit the given number of lipids in.
        """
        parser.add_argument("-pbc", required=False, type=str, default="hexagonal", help="PBC type: hexagonal, rectangular, square, cubic, optimal or keep")
        parser.add_argument("-d", required=False, type=float, default=0, help="Distance between periodic images (nm)")
        parser.add_argument("-dz", required=False, type=float, default=0, help="Z distance between periodic images (nm)")
        parser.add_argument("-x", required=False, type=vector, default=0, help="X dimension or first lattice vector of system (nm)")
        parser.add_argument("-y", required=False, type=vector, default=0, help="Y dimension or first lattice vector of system (nm)")
        parser.add_argument("-z", required=False, type=vector, default=0, help="Z dimension or first lattice vector of system (nm)")
        parser.add_argument("-box", required=False, type=read_box, default=None, help="Box in GRO (3 or 9 floats) or PDB (6 floats) format, comma separated")
        parser.add_argument("-n", required=False, type=str, default=None, help="Index file --- TO BE IMPLEMENTED")
        # Membrane / Lipid options
        """
            Membrane/lipid related options.
            The options -l and -u can be given multiple times. Option -u can be
            used to set the lipid type and abundance for the upper leaflet. Option
            -l sets the type and abundance for the lower leaflet if option -u is
            also given, or for both leaflets if option -u is not given. The
            meaning of the number depends on whether option -d is used to set up
            PBC
        """
        parser.add_argument("-l", required=False, nargs="+", type=str, default=None, action="append", help="Lipid type and relative abundance in lower leaf (NAME[:#])")
        parser.add_argument("-u", required=False, nargs="+", type=str, default=None, action="append", help="Lipid type and relative abundance in upper leaf (NAME[:#])")
        parser.add_argument("-a", required=False, type=float, default=0.60, help="Area per lipid (nm*nm)")
        parser.add_argument("-au", required=False, type=float, default=None, help="Area per lipid (nm*nm) for upper layer")
        parser.add_argument("-asym", required=False, type=int, default=None, help="Membrane asymmetry (number of lipids)")
        parser.add_argument("-hole", required=False, type=float, default=None, help="Make a hole in the membrane with specified radius")
        parser.add_argument("-disc", required=False, type=float, default=None, help="Make a membrane disc with specified radius")
        parser.add_argument("-rand", required=False, type=float, default=0.1, help="Random kick size (maximum atom displacement)")
        parser.add_argument("-bd", required=False, type=float, default=0.3, help="Bead distance unit for scaling z-coordinates (nm)")
        # Protein related options
        parser.add_argument("-center", required=False, default=False, action="store_true", help="Center the protein on z")
        parser.add_argument("-orient", required=False, default=False, action="store_true", help="Orient protein in membrane")
        parser.add_argument("-rotate", required=False, type=str, default=None, help="Rotate protein (random|princ|angle(float))")
        parser.add_argument("-od", required=False, type=float, default=1.0, help="Grid spacing for determining orientation")
        parser.add_argument("-op", required=False, type=float, default=4.0, help="Hydrophobic ratio power for determining orientation")
        parser.add_argument("-fudge", required=False, type=float, default=0.1, help="Fudge factor for allowing lipid-protein overlap")
        parser.add_argument("-ring", required=False, default=False, action="store_true", help="Put lipids inside the protein")
        parser.add_argument("-dm", required=False, type=float, default=None, help="Shift protein with respect to membrane")
        # Solvent related options
        parser.add_argument("-sol", required=False, nargs="+", type=str, default=None, action="append", help="Solvent type and relative abundance (NAME[:#])")
        parser.add_argument("-sold", required=False, type=float, default=0.5, help="Solvent diameter")
        parser.add_argument("-solr", required=False, type=float, default=0.1, help="Solvent random kick")
        parser.add_argument("-excl", required=False, type=float, default=1.5, help="Exclusion range (nm) for solvent addition relative to membrane center")
        # Salt related options
        parser.add_argument("-salt", required=False, type=str, default=None, help="Salt concentration")
        parser.add_argument("-charge", required=False, type=str, default="auto", help="Charge of system. Set to auto to infer from residue names")
        # Define additional lipid types (same format as in lipid-martini-itp-v01.py)
        parser.add_argument("-alname", required=False, type=str, nargs="+", default=None, action="append", help="Additional lipid name, x4 letter")
        parser.add_argument("-alhead", required=False, type=str, nargs="+", default=None, action="append", help="Additional lipid head specification string")
        parser.add_argument("-allink", required=False, type=str, nargs="+", default=None, action="append", help="Additional lipid linker specification string")
        parser.add_argument("-altail", required=False, type=str, nargs="+", default=None, action="append", help="Additional lipid tail specification string")
        parser.add_argument("-alcharge", required=False, type=str, nargs="+", default=None, action="append", help="Additional lipid charge")
        parser.add_argument("-m", required=False, type=str, nargs="+", default=None, action="append", help="Read molecule definitions from file(s)")
        # Parse arguments
        args = parser.parse_args()
        # Translate args
        self.options = vars(args)
        # Unpack multi arguments from [[x], ...] to [x, ...]
        multi_args = ("alname", "alhead", "allink", "altail", "alcharge", "m", "sol", "l", "u", "f")
        for arg in multi_args:
            self.options[arg] = [_[0] for _ in self.options[arg]] if self.options[arg] is not None else []
        # Map arguments to their explicit variables
        self.tm = self.options["f"]
        self.lipL = self.options["l"]
        self.lipU = self.options["u"]
        self.solv = self.options["sol"]
        self.usrnames = self.options["alname"]
        self.usrheads = self.options["alhead"]
        self.usrlinks = self.options["allink"]
        self.usrtails = self.options["altail"]
        self.usrcharg = self.options["alcharge"]
        self.mollist = self.options["m"]
        print(self.options)

    def process_solvents(self):
        # Charge of the system so far

        last = None
        mcharge = 0
        for j in self.membrane.atoms:
            if not j[0].strip().startswith('v') and j[1:3] != last:
                mcharge += charges.get(j[1].strip(), 0)
            last = j[1:3]

        last = None
        pcharge = 0
        for j in self.protein.atoms:
            if not j[0].strip().startswith('v') and j[1:3] != last:
                pcharge += charges.get(j[1].strip(), 0)
            last = j[1:3]

        # mcharge = sum([charges.get(i[0].strip(),0) for i in set([j[1:3] for j in membrane.atoms])])
        # pcharge = sum([charges.get(i[0].strip(),0) for i in set([j[1:3] for j in protein.atoms if not j[0].strip().startswith('v')])])

        charge = mcharge + pcharge
        plen, mlen, slen = 0, 0, 0
        plen = self.protein and len(self.protein) or 0
        print("; NDX Solute {} {}".format(1, self.protein and plen or 0), file=sys.stderr)
        print("; Charge of protein: {}".format(pcharge), file=sys.stderr)

        mlen = self.membrane and len(self.membrane) or 0
        print("; NDX Membrane {} {}".format(1 + plen, self.membrane and plen + mlen or 0), file=sys.stderr)
        print("; Charge of membrane: {}".format(mcharge), file=sys.stderr)
        print("; Total charge: {}".format(charge), file=sys.stderr)

        if self.solv:

            # Set up a grid
            d = 1 / self.options["sold"]

            nx, ny, nz = int(1 + d * self.pbcx), int(1 + d * self.pbcy), int(1 + d * self.pbcz)
            dx, dy, dz = self.pbcx / nx, self.pbcy / ny, self.pbcz / nz
            excl, hz = int(nz * self.options["excl"] / self.pbcz), int(0.5 * nz)

            zshift = 0
            if self.membrane:
                memz = [i[2] for i in self.membrane.coord]
                midz = (max(memz) + min(memz)) / 2
                hz = int(nz * midz / self.pbcz)  # Grid layer in which the membrane is located
                zshift = (hz + 0.5) * nz - midz  # Shift of membrane middle to center of grid layer

            # Initialize a grid of solvent, spanning the whole cell
            # Exclude all cells within specified distance from membrane center
            grid = [[[i < hz - excl or i > hz + excl for i in range(nz)] for j in range(ny)] for i in range(nx)]

            # Flag all cells occupied by protein or membrane
            for p, q, r in self.protein.coord + self.membrane.coord:
                for _s, t, u in Insane3.points_on_sphere(20):
                    x, y, z = p + 0.33 * _s, q + 0.33 * t, r + 0.33 * u
                    if z >= self.pbcz:
                        x -= self.box[2][0]
                        y -= self.box[2][1]
                        z -= self.box[2][2]
                    if z < 0:
                        x += self.box[2][0]
                        y += self.box[2][1]
                        z += self.box[2][2]
                    if y >= self.pbcy:
                        x -= self.box[1][0]
                        y -= self.box[1][1]
                    if y < 0:
                        x += self.box[1][0]
                        y += self.box[1][1]
                    if x >= self.pbcx:
                        x -= self.box[0][0]
                    if x < 0:
                        x += self.box[0][0]
                    grid[int(nx * x / self.rx)][int(ny * y / self.ry)][int(nz * z / self.rz)] = False

            # Set the center for each solvent molecule
            kick = self.options["solr"]
            grid = [(R(), (i + 0.5 + R() * kick) * dx, (j + 0.5 + R() * kick) * dy, (k + 0.5 + R() * kick) * dz)
                    for i in range(nx) for j in range(ny) for k in range(nz) if grid[i][j][k]]

            # Sort on the random number
            grid.sort()

            # 'grid' contains all positions on which a solvent molecule can be placed.
            # The number of positions is taken as the basis for determining the salt concentration.
            # This is fine for simple salt solutions, but may not be optimal for complex mixtures
            # (like when mixing a 1M solution of this with a 1M solution of that

            # First get names and relative numbers for each solvent
            solnames, solnums = list(zip(*[parse_mol(i) for i in self.solv]))
            solnames, solnums = list(solnames), list(solnums)
            totS = float(sum(solnums))

            # Set the number of ions to add
            nna, ncl = 0, 0
            if self.options["salt"]:

                # If the concentration is set negative, set the charge to zero
                if self.options["salt"].startswith("-"):
                    charge = 0
                    self.options["salt"] = -float(self.options["salt"])
                else:
                    self.options["salt"] = float(self.options["salt"])

                # Determine charge to use, either determined or given on command line
                if self.options["charge"] != "0":
                    charge = (self.options["charge"] != "auto") and int(self.options["charge"]) or charge
                else:
                    charge = 0

                # Determine number of sodium and chloride to add
                concentration = self.options["salt"]
                nsol = ("SPC" in solnames and 1 or 4) * len(grid)
                ncl = max(max(0, charge), int(.5 + .5 * (concentration * nsol / (27.7 + concentration) + charge)))
                nna = ncl - charge

            # Correct number of grid cells for placement of solvent
            ngrid = len(grid) - nna - ncl
            num_sol = [int(ngrid * i / totS) for i in solnums]

            # Add salt to solnames and num_sol
            if nna:
                solnames.append("NA+")
                num_sol.append(nna)
                self.solv.append("NA+")
            if ncl:
                solnames.append("CL-")
                num_sol.append(ncl)
                self.solv.append("CL-")

            # Names and grid positions for solvent molecules
            solvent = list(zip([_s for i, _s in list(zip(num_sol, solnames)) for j in range(i)], grid))

            # Extend the list of molecules (for the topology)
            self.molecules.extend(list(zip(solnames, num_sol)))

            # Build the solvent
            for resn, (rndm, x, y, z) in solvent:
                self.resi += 1
                solmol = solventParticles.get(resn)
                if solmol and len(solmol) > 1:
                    # Random rotation (quaternion)
                    u, v, w = random.random(), 2 * math.pi * random.random(), 2 * math.pi * random.random()
                    s, t = math.sqrt(1 - u), math.sqrt(u)
                    qw, qx, qy, qz = s * math.sin(v), s * math.cos(v), t * math.sin(w), t * math.cos(w)
                    qq = qw * qw - qx * qx - qy * qy - qz * qz
                    for atnm, (px, py, pz) in solmol:
                        qp = 2 * (qx * px + qy * py + qz * pz)
                        rx = x + qp * qx + qq * px + qw * (qy * pz - qz * py)
                        ry = y + qp * qy + qq * py + qw * (qz * px - qx * pz)
                        rz = z + qp * qz + qq * pz + qw * (qx * py - qy * px)
                        self.sol.append(("%5d%-5s%5s%5d" % (self.resi % 1e5, resn, atnm, self.atid % 1e5), (rx, ry, rz)))
                        self.atid += 1
                else:
                    self.sol.append(
                        ("%5d%-5s%5s%5d" % (self.resi % 1e5, resn, solmol and solmol[0][0] or resn, self.atid % 1e5), (x, y, z)))
                    self.atid += 1

        else:
            solvent, sol = None, []

        slen = solvent and len(self.sol) or 0
        print("; NDX Solvent {} {}".format(1 + plen + mlen, solvent and plen + mlen + slen or 0), file=sys.stderr)
        print("; NDX System {} {}".format(1, plen + mlen + slen), file=sys.stderr)
        print("; \"I mean, the good stuff is just INSANE\" --Julia Ormond", file=sys.stderr)

    def process_membrane(self):
        """
        Spawn a membrane.
        :return: None
        """
        if self.lipL:
            # Lipids are added on grid positions, using the prototypes defined above.
            # If a grid position is already occupied by protein, the position is untagged.

            lipd = self.lo_lipd

            # Number of lipids in x and y in lower leaflet if there were no solute
            lo_lipids_x = int(self.pbcx / lipd + 0.5)
            lo_lipdx = self.pbcx / lo_lipids_x
            lo_rlipx = range(lo_lipids_x)
            lo_lipids_y = int(self.pbcy / lipd + 0.5)
            lo_lipdy = self.pbcy / lo_lipids_y
            lo_rlipy = range(lo_lipids_y)

            if self.options["au"]:
                lipd = self.up_lipd

            # Number of lipids in x and y in upper leaflet if there were no solute
            up_lipids_x = int(self.pbcx / lipd + 0.5)
            up_lipdx = self.pbcx / up_lipids_x
            up_rlipx = range(up_lipids_x)
            up_lipids_y = int(self.pbcy / lipd + 0.5)
            up_lipdy = self.pbcy / up_lipids_y
            up_rlipy = range(up_lipids_y)

            # Set up grids to check where to place the lipids
            grid_lo = [[0 for j in lo_rlipy] for i in lo_rlipx]
            grid_up = [[0 for j in up_rlipy] for i in up_rlipx]

            # If there is a protein, mark the corresponding cells as occupied
            if self.protein:
                # Calculate number density per cell
                for i in self.prot_lo:
                    grid_lo[int(lo_lipids_x * i[0] / self.rx) % lo_lipids_x][int(lo_lipids_y * i[1] / self.ry) % lo_lipids_y] += 1
                for i in self.prot_up:
                    grid_up[int(up_lipids_x * i[0] / self.rx) % up_lipids_x][int(up_lipids_y * i[1] / self.ry) % up_lipids_y] += 1

            # Determine which cells to consider occupied, given the fudge factor
            # The array is changed to boolean type here
            maxd = float(max([max(i) for i in grid_up + grid_lo]))
            if maxd == 0:
                if self.protein:
                    print("; The protein seems not to be inside the membrane.", file=sys.stderr)
                    print("; Run with -orient to put it in.", file=sys.stderr)
                maxd = 1

            fudge = self.options["fudge"]
            grid_up = [[(j / maxd) <= fudge for j in i] for i in grid_up]
            grid_lo = [[(j / maxd) <= fudge for j in i] for i in grid_lo]

            # If we don't want lipids inside of the protein
            # we also mark everything from the center up to the first cell filled
            if not self.options["ring"]:

                # Upper leaflet
                marked = [(i, j) for i in up_rlipx for j in up_rlipy if not grid_up[i][j]]
                if marked:
                    # Find the center
                    cx, cy = [float(sum(i)) / len(marked) for i in list(zip(*marked))]
                    for i, j in marked:
                        md = int(abs(i - cx) + abs(j - cy))  # Manhattan length/distance
                        for f in range(md):
                            ii = int(cx + f * (i - cx) / md)
                            jj = int(cy + f * (j - cy) / md)
                            grid_up[ii][jj] = False

                # Lower leaflet
                marked = [(i, j) for i in lo_rlipx for j in lo_rlipy if not grid_lo[i][j]]
                if marked:
                    # Find the center
                    cx, cy = [float(sum(i)) / len(marked) for i in list(zip(*marked))]
                    for i, j in marked:
                        md = int(abs(i - cx) + abs(j - cy))  # Manhattan length
                        for f in range(md):
                            ii = int(cx + f * (i - cx) / md)
                            jj = int(cy + f * (j - cy) / md)
                            grid_lo[ii][jj] = False

                # If we make a circular patch, we flag the cells further from the
                # protein or box center than the given radius as occupied.
                if self.options["disc"]:
                    if self.protein:
                        cx, cy = self.protein.center()[:2]
                    else:
                        cx, cy = 0.5 * self.pbcx, 0.5 * self.pbcy
                    for i in range(len(grid_lo)):
                        for j in range(len(grid_lo[i])):
                            if (i * self.pbcx / lo_lipids_x - cx) ** 2 + (j * self.pbcy / lo_lipids_y - cy) ** 2 > self.options["disc"] ** 2:
                                grid_lo[i][j] = False
                    for i in range(len(grid_up)):
                        for j in range(len(grid_up[i])):
                            if (i * self.pbcx / up_lipids_x - cx) ** 2 + (j * self.pbcy / up_lipids_y - cy) ** 2 > self.options["disc"] ** 2:
                                grid_up[i][j] = False

                # If we need to add a hole, we simply flag the corresponding cells
                # as occupied. The position of the hole depends on the type of PBC,
                # to ensure an optimal arrangement of holes around the protein. If
                # there is no protein, the hole is just put in the center.
                if self.options["hole"]:
                    # Lower leaflet
                    if self.protein:
                        if ("square".startswith(self.options["pbc"]) or
                                "rectangular".startswith(self.options["pbc"])):
                            hx, hy = (0, 0)
                        else:
                            hx, hy = (0, int(lo_lipids_y * math.cos(math.pi / 6) / 9 + 0.5))
                    else:
                        hx, hy = (int(0.5 * lo_lipids_x), int(0.5 * lo_lipids_y))
                    hr = int(self.options["hole"] / min(lo_lipdx, lo_lipdy) + 0.5)
                    ys = int(lo_lipids_x * self.box[1][0] / self.box[0][0] + 0.5)
                    print("; Making a hole with radius {} nm centered at grid cell ({},{})".format(self.options["hole"], hx, hy), hr, file=sys.stderr)
                    hr -= 1
                    for ii in range(hx - hr - 1, hx + hr + 1):
                        for jj in range(hx - hr - 1, hx + hr + 1):
                            xi, yj = ii, jj
                            if (ii - hx) ** 2 + (jj - hy) ** 2 < hr ** 2:
                                if jj < 0:
                                    xi += ys
                                    yj += lo_lipids_y
                                if jj >= lo_lipids_y:
                                    xi -= ys
                                    yj -= lo_lipids_y
                                if xi < 0:
                                    xi += lo_lipids_x
                                if xi >= lo_lipids_x:
                                    xi -= lo_lipids_x
                                grid_lo[xi][yj] = False
                                grid_up[xi][yj] = False
                    # Upper leaflet
                    if self.protein:
                        if ("square".startswith(self.options["pbc"]) or
                                "rectangular".startswith(self.options["pbc"])):
                            hx, hy = (0, 0)
                        else:
                            hx, hy = (0, int(up_lipids_y * math.cos(math.pi / 6) / 9 + 0.5))
                    else:
                        hx, hy = (int(0.5 * up_lipids_x), int(0.5 * up_lipids_y))
                    hr = int(self.options["hole"] / min(up_lipdx, up_lipdy) + 0.5)
                    ys = int(up_lipids_x * self.box[1][0] / self.box[0][0] + 0.5)
                    print("; Making a hole with radius {} nm centered at grid cell ({},{})".format(self.options["hole"], hx, hy), hr, file=sys.stderr)
                    hr -= 1
                    for ii in range(hx - hr - 1, hx + hr + 1):
                        for jj in range(hx - hr - 1, hx + hr + 1):
                            xi, yj = ii, jj
                            if (ii - hx) ** 2 + (jj - hy) ** 2 < hr ** 2:
                                if jj < 0:
                                    xi += ys
                                    yj += up_lipids_y
                                if jj >= up_lipids_y:
                                    xi -= ys
                                    yj -= up_lipids_y
                                if xi < 0:
                                    xi += up_lipids_x
                                if xi >= up_lipids_x:
                                    xi -= up_lipids_x
                                grid_up[xi][yj] = False

                # Set the XY coordinates
                # To randomize the lipids we add a random number which is used for sorting
                upper, lower = [], []
                for i in range(up_lipids_x):
                    for j in range(up_lipids_y):
                        if grid_up[i][j]:
                            upper.append((random.random(), i * self.pbcx / up_lipids_x, j * self.pbcy / up_lipids_y))
                for i in range(lo_lipids_x):
                    for j in range(lo_lipids_y):
                        if grid_lo[i][j]:
                            lower.append((random.random(), i * self.pbcx / lo_lipids_x, j * self.pbcy / lo_lipids_y))

                # Sort on the random number
                upper.sort()
                lower.sort()

                # Extract coordinates, taking asymmetry in account
                asym = self.options["asym"] or 0
                upper = [i[1:] for i in upper[max(0, asym):]]
                lower = [i[1:] for i in lower[max(0, -asym):]]

                print("; X: {:5.3f} ({} bins) Y: {:5.3f} ({} bins) in upper leaflet".format(self.pbcx, up_lipids_x, self.pbcy, up_lipids_y), file=sys.stderr)
                print("; X: {:5.3f} ({} bins) Y: {:5.3f} ({} bins) in lower leaflet".format(self.pbcx, lo_lipids_x, self.pbcy, lo_lipids_y), file=sys.stderr)
                print("; {} lipids in upper leaflet, {} lipids in lower leaflet".format(len(upper), len(lower)), file=sys.stderr)

                # Types of lipids, relative numbers, fractions and numbers

                self.lipU = self.lipU or self.lipL

                # Upper leaflet (+1)
                self.lipU, self.numU = list(zip(*[parse_mol(i) for i in self.lipU]))
                totU = float(sum(self.numU))
                num_up = [int(len(upper) * i / totU) for i in self.numU]
                lip_up = [l for i, l in list(zip(num_up, self.lipU)) for j in range(i)]
                leaf_up = (1, list(zip(lip_up, upper)), self.up_lipd, up_lipdx, up_lipdy)

                # Lower leaflet (-1)
                self.lipL, self.numL = list(zip(*[parse_mol(i) for i in self.lipL]))
                totL = float(sum(self.numL))
                num_lo = [int(len(lower) * i / totL) for i in self.numL]
                lip_lo = [l for i, l in list(zip(num_lo, self.lipL)) for j in range(i)]
                leaf_lo = (-1, list(zip(lip_lo, lower)), self.lo_lipd, lo_lipdx, lo_lipdy)

                self.molecules = list(zip(self.lipU, num_up)) + list(zip(self.lipL, num_lo))

                kick = self.options["rand"]

                # Build the membrane
                for leaflet, leaf_lip, lipd, lipdx, lipdy in [leaf_up, leaf_lo]:
                    for lipid, pos in leaf_lip:
                        # Increase the residue number by one
                        self.resi += 1
                        # Set the random rotation for this lipid
                        rangle = 2 * random.random() * math.pi
                        rcos = math.cos(rangle)
                        rsin = math.sin(rangle)
                        rcosx = rcos * lipdx * 2 / 3
                        rcosy = rcos * lipdy * 2 / 3
                        rsinx = rsin * lipdx * 2 / 3
                        rsiny = rsin * lipdy * 2 / 3
                        # Fetch the atom list with x,y,z coordinates
                        # atoms    = zip(lipidsa[lipid][1].split(),lipidsx[lipidsa[lipid][0]],lipidsy[lipidsa[lipid][0]],lipidsz[lipidsa[lipid][0]])
                        # Only keep atoms appropriate for the lipid
                        # at,ax,ay,az = zip(*[i for i in atoms if i[0] != "-"])
                        at, ax, ay, az = list(zip(*self.liplist[lipid].build(diam=lipd)))
                        # The z-coordinates are spaced at 0.3 nm,
                        # starting with the first bead at 0.15 nm
                        az = [leaflet * (0.5 + (i - min(az))) * self.options["bd"] for i in az]
                        xx = list(zip(ax, ay))
                        nx = [rcosx * i - rsiny * j + pos[0] + lipdx / 2 + random.random() * kick for i, j in xx]
                        ny = [rsinx * i + rcosy * j + pos[1] + lipdy / 2 + random.random() * kick for i, j in xx]
                        # Add the atoms to the list
                        for i in range(len(at)):
                            atom = "%5d%-5s%5s%5d" % (self.resi, lipid, at[i], self.atid)
                            self.membrane.coord.append((nx[i], ny[i], az[i]))
                            self.membrane.atoms.append((at[i], lipid, self.resi, 0, 0, 0))
                            self.atid += 1

                # Now move everything to the center of the box before adding solvent
                mz = self.pbcz / 2
                z = [i[2] for i in self.protein.coord + self.membrane.coord]
                mz -= (max(z) + min(z)) / 2
                self.protein += (0, 0, mz)
                self.membrane += (0, 0, mz)

    def process_lipids(self):
        """
        Add lipids.
        :return: Array<Lipid>
        """
        liplist = LipidList()
        # First add internal lipids
        for name, lip in lipidsa.items():
            _moltype = lip[0]
            template = list(zip(lipidsx[_moltype], lipidsy[_moltype], lipidsz[_moltype]))
            liplist[name] = Lipid(name=name, beads=lip[1], template=template)

        # Then add lipids from file
        for filename in self.mollist:
            stuff = ""
            with open(filename) as _fh:
                stuff = _fh.read().split("@INSANE")
            for group in stuff[1:]:
                lines = group.split("\n")
                lipdef = lines.pop(0)
                beads = None
                for line in lines:
                    if line.startswith('[') or not line.strip():
                        break
                    if "@BEADS" in line:
                        beads = line.split("@BEADS")[1].split()
                lip = Lipid(string=lipdef, beads=beads)
                liplist[lip.name] = lip

        # Last, add lipids from command line
        for name, head, link, tail in list(zip(self.usrnames, self.usrheads, self.usrlinks, self.usrtails)):
            heads = head.replace(".", " ").split()
            linkers = link.replace(".", " ").split()
            tails = tail.replace(".", " ").split()
            liplist[name] = Lipid(name=name, head=heads, link=linkers, tail=tails)
        return liplist

    def process_box(self):
        """
        Set periodic boundary conditions.
        :return: None
        """
        # option -box overrides everything
        if self.options["box"]:
            self.options["x"] = self.options["box"][:3]
            self.options["y"] = self.options["box"][3:6]
            self.options["z"] = self.options["box"][6:]
        # option -pbc keep really overrides everything
        if self.options["pbc"] == "keep" and self.tm:
            self.options["x"] = self.tm[0].box[:3]
            self.options["y"] = self.tm[0].box[3:6]
            self.options["z"] = self.tm[0].box[6:]
        # options -x, -y, -z take precedence over automatic determination
        if isinstance(self.options["x"], (list, tuple)):
            self.pbc_set_x = self.options["x"]
        elif self.options["x"]:
            self.pbc_set_x = [self.options["x"], 0, 0]

        if isinstance(self.options["y"], (list, tuple)):
            self.pbc_set_y = self.options["y"]
        elif self.options["y"]:
            self.pbc_set_y = [0, self.options["y"], 0]

        if isinstance(self.options["z"], (list, tuple)):
            self.pbc_set_z = self.options["z"]
        elif self.options["z"]:
            self.pbc_set_z = [0, 0, self.options["z"]]

    def process_protein(self):
        """
        Add and orient protein if available.
        :return: None
        """
        if not self.tm:
            self.__process_membrane_no_protein()
        else:
            for prot in self.tm:
                if not self.lipL:
                    prot = self.__process_protein_no_membrane(prot)
                else:
                    prot = self.__process_protein_with_membrane(prot)
                # And we collect the atoms
                self.protein.atoms.extend(prot.atoms)
                self.protein.coord.extend(prot.coord)

            # Extract the parts of the protein that are in either leaflet
            for ix, iy, iz in self.protein.coord:
                if 0 < iz < 2.4:
                    self.prot_up.append((ix, iy))
                elif 0 > iz > -2.4:
                    self.prot_lo.append((ix, iy))

            # Current residue ID is set to that of the last atom
            self.resi = self.protein.atoms[-1][2]

        self.atid = len(self.protein) + 1

        # The box dimensions are now (likely) set.
        # If a protein was given, it is positioned in the center of the
        # rectangular brick.

        # Set the lattice vectors
        if ("rectangular".startswith(self.options["pbc"]) or
                "square".startswith(self.options["pbc"]) or
                "cubic".startswith(self.options["pbc"])):
            self.box = [[self.pbcx, 0, 0], [0, self.pbcy, 0], [0, 0, self.pbcz]]
        elif not self.lipL:
            # Rhombic dodecahedron with square XY plane
            self.box = [[self.pbcx, 0, 0], [0, self.pbcy, 0], [0.5 * self.pbcx, 0.5 * self.pbcx, self.pbcz]]
        elif "hexagonal".startswith(self.options["pbc"]):
            self.box = [[self.pbcx, 0, 0], [math.sin(math.pi / 6) * self.pbcx, self.pbcy, 0], [0, 0, self.pbcz]]
        else:  # optimal packing; rhombic dodecahedron with hexagonal XY plane
            self.box = [[self.pbcx, 0, 0], [math.sin(math.pi / 6) * self.pbcx, self.pbcy, 0], [self.pbcx / 2, self.pbcy / 3, self.pbcz]]

        # Override lattice vectors if they were set explicitly
        self.box[0] = self.pbc_set_x or self.box[0]
        self.box[1] = self.pbc_set_y or self.box[1]
        self.box[2] = self.pbc_set_z or self.box[2]

        self.grobox = (
            self.box[0][0], self.box[1][1], self.box[2][2],
            self.box[0][1], self.box[0][2], self.box[1][0],
            self.box[1][2], self.box[2][0], self.box[2][1]
        )

        self.pbcx, self.pbcy, self.pbcz = self.box[0][0], self.box[1][1], self.box[2][2]

        self.rx, self.ry, self.rz = self.pbcx + 1e-8, self.pbcy + 1e-8, self.pbcz + 1e-8

    def __process_protein_with_membrane(self, prot):
        """
        Build membrane around a protein.
        :param prot: Structure
        :return: Structure
        """
        # Have to build a membrane around the protein.
        # So first put the protein in properly.

        # Center the protein and store the shift
        shift = prot.center((0, 0, 0))

        # 1. Orient with respect to membrane
        # Orient the protein according to the TM region, if requested
        # This doesn't actually work very well...
        if self.options["orient"]:

            # Grid spacing (nm)
            d = self.options["od"]
            pw = self.options["op"]

            # Determine grid size
            mx, my, mz = prot.fun(min)
            rx, ry, rz = prot.fun(lambda x: max(x) - min(x) + 1e-8)

            # Number of grid cells
            nx, ny, nz = int(rx / d + 0.5), int(ry / d + 0.5), int(rz / d + 0.5)

            # Initialize grids
            atom = [[[0 for i in range(nz + 2)] for j in range(ny + 2)] for k in range(nx + 2)]
            phobic = [[[0 for i in range(nz + 2)] for j in range(ny + 2)] for k in range(nx + 2)]
            surface = []
            for i, (ix, iy, iz) in zip(prot.atoms, prot.coord):
                if i[1] != "DUM":
                    jx, jy, jz = int(nx * (ix - mx) / rx), int(ny * (iy - my) / ry), int(nz * (iz - mz) / rz)
                    atom[jx][jy][jz] += 1
                    phobic[jx][jy][jz] += (i[1].strip() in apolar)

            # Determine average density
            occupd = sum([bool(k) for i in atom for j in i for k in j])
            avdens = float(sum([sum(j) for i in atom for j in i])) / occupd

            # cgofile  = open('density.cgo',"w")
            # cgofile.write('[\n')
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if atom[i][j][k] > 0.1 * avdens:
                            # Check the neighbouring cells; If one of them is not occupied, count cell as surface
                            if not (atom[i - 1][j][k] and atom[i + 1][j][k] and
                                        atom[i][j - 1][k] and atom[i][j + 1][k] and
                                        atom[i][j][k - 1] and atom[i][j][k + 1]):
                                sx, sy, sz = mx + rx * (i + 0.5) / nx, my + ry * (j + 0.5) / ny, mz + rz * (
                                k + 0.5) / nz
                                sw = (2.0 * phobic[i][j][k] / atom[i][j][k]) ** pw
                                surface.append((sx, sy, sz, sw))
                                # cgofile.write("    7.0, %f, %f, %f, %f,\n"%(10*sx,10*sy,10*sz,0.25*sw))
            # cgofile.write(']\n')
            # cgofile.close()

            sx, sy, sz, w = zip(*surface)
            big_w = 1.0 / sum(w)

            # Weighted center of apolar region; has to go to (0,0,0)
            sxm, sym, szm = [sum(p) * big_w for p in zip(*[(m * i, m * j, m * k) for m, i, j, k in zip(w, sx, sy, sz)])]

            # Place apolar center at origin
            prot.center((-sxm, -sym, -szm))
            sx, sy, sz = zip(*[(i - sxm, j - sym, k - szm) for i, j, k in zip(sx, sy, sz)])

            # Determine weighted deviations from centers
            dx, dy, dz = zip(*[(m * i, m * j, m * k) for m, i, j, k in zip(w, sx, sy, sz)])

            # Covariance matrix for surface
            xx, yy, zz, xy, yz, zx = [sum(p) * big_w for p in
                                      zip(*[(i * i, j * j, k * k, i * j, j * k, k * i) for i, j, k in zip(dx, dy, dz)])]

            # PCA: u,v,w are a rotation matrix
            (ux, uy, uz), (vx, vy, vz), (wx, wy, wz), r = mijn_eigen_sym_3x3(xx, yy, zz, xy, zx, yz)

            # Rotate the coordinates
            prot.coord = [(ux * i + uy * j + uz * k, vx * i + vy * j + vz * k, wx * i + wy * j + wz * k) for i, j, k in
                          prot.coord]

        # 4. Orient the protein in the xy-plane
        # i. According to principal axes and unit cell
        if self.options["rotate"] == "princ":

            x, y, z = zip(*prot.coord)

            # The rotation matrix in the plane equals the transpose
            # of the matrix of eigenvectors from the 2x2 covariance
            # matrix of the positions.
            # For numerical stability we do
            # d_i     = x_i - x_0
            # mean(x) = x_0 + sum(d_i)/N =
            # var(x)  = sum((d_i - mean(d))**2)/(N-1)
            xy = ssd(x, y)
            if xy != 0:
                xx = ssd(x, x)
                yy = ssd(y, y)

                # The eigenvalues are the roots of the 2nd order
                # characteristic polynomial, with the coefficients
                # equal to the trace and the determinant of the
                # matrix.
                t, d = xx + yy, xx * yy - xy * xy
                # The two eigenvectors form a 2D rotation matrix
                # R = ((cos,sin),(-sin,cos)), which means that
                # the second eigenvector follows directly from
                # the first. We thus only need to determine one.
                l1 = t / 2 + math.sqrt(0.25 * t * t - d)

                ux, uy = l1 - yy, xy
                lu = math.sqrt(ux * ux + uy * uy)

                ux /= lu
                uy /= lu

                # Finally we rotate the system in the plane by
                # matrix multiplication with the transpose of
                # the matrix of eigenvectors
                prot.coord = [(ux * i + uy * j, ux * j - uy * i, k) for i, j, k in zip(x, y, z)]

        # ii. Randomly
        elif self.options["rotate"] == "random":
            ux = math.cos(R() * 2 * math.pi)
            uy = math.sqrt(1 - ux * ux)
            prot.coord = [(ux * i + uy * j, ux * j - uy * i, k) for i, j, k in prot.coord]

        # iii. Specifically
        elif self.options["rotate"]:
            ux = math.cos(float(self.options["rotate"]) * math.pi / 180.)
            uy = math.sin(float(self.options["rotate"]) * math.pi / 180.)
            prot.coord = [(ux * i + uy * j, ux * j - uy * i, k) for i, j, k in prot.coord]

        # 5. Determine the minimum and maximum x and y of the protein
        pmin, pmax = prot.fun(min), prot.fun(max)
        prng = (pmax[0] - pmin[0], pmax[1] - pmin[1], pmax[2] - pmin[2])
        center = (0.5 * (pmin[0] + pmax[0]), 0.5 * (pmin[1] + pmax[1]))

        # Set the z-dimension
        self.pbcz = self.pbc_set_z and self.pbc_set_z[2]
        # If it is not set, set pbcz to the dimension of the protein
        self.pbcz = self.pbcz or prng[2]
        self.pbcz += self.options["dz"] or self.options["d"] or 0

        # At this point we should shift the subsequent proteins such that they end up
        # at the specified distance, in case we have a number of them to do
        # y-shift is always -ycenter
        # x-shift is -xmin+distance+xmax(current)
        xshft, yshft = self.xshifts[-1] - pmin[0] + (self.options["d"] or 0), -center[1]
        self.xshifts.append(self.xshifts[-1] + pmax[0] + (self.options["d"] or 0))

        # 6. Set box (brick) dimensions
        if self.options["disc"]:
            self.pbcx = self.options["d"] + 2 * self.options["disc"]
            if ("square".startswith(self.options["pbc"]) or
                    "rectangular".startswith(self.options["pbc"])):
                self.pbcy = self.pbcx
            else:
                self.pbcy = math.cos(math.pi / 6) * self.pbcx
        else:
            self.pbcx = (self.options["d"] or 0) + prng[0]
            if "square".startswith(self.options["pbc"]):
                self.pbcy = self.pbcx
            elif "rectangular".startswith(self.options["pbc"]):
                self.pbcy = self.options["d"] + prng[1]
            else:
                # This goes for a hexagonal cell as well as for the optimal arrangement
                # The latter is hexagonal in the membrane plane anyway...
                self.pbcy = math.cos(math.pi / 6) * self.pbcx

        # 7. Adjust PBC for hole
        # If we need to add a hole, we have to scale the system
        # The scaling depends on the type of PBC
        if self.options["hole"]:
            if ("square".startswith(self.options["pbc"]) or
                    "rectangular".startswith(self.options["pbc"])):
                scale = 1 + self.options["hole"] / min(self.pbcx, self.pbcy)
            else:
                area = self.options["hole"] ** 2 / math.cos(math.pi / 6)
                scale = 1 + area / (self.pbcx * self.pbcy)
            self.pbcx, self.pbcy = scale * self.pbcx, scale * self.pbcy

        self.pbcx = self.pbc_set_x and self.pbc_set_x[0] or self.pbcx
        self.pbcy = self.pbc_set_y and self.pbc_set_y[1] or self.pbcy

        # 2. Shift of protein relative to the membrane center
        zshift = 0
        if not self.options["center"]:
            zshift = -shift[2]
        if self.options["dm"]:
            if self.options["dm"] < 0:
                zshift += self.options["dm"]  # - max(zip(*prot.coord)[2])
            else:
                zshift += self.options["dm"]  # - min(zip(*prot.coord)[2])

        # Now we center the system in the rectangular
        # brick corresponding to the unit cell
        # If -center is given, also center z in plane
        prot += (0.5 * self.pbcx, 0.5 * self.pbcy, zshift)
        return prot

    def __process_protein_no_membrane(self, prot):
        """
        Solvate protein, do not add a membrane.
        :param prot: Structure
        :return: Structure
        """
        # A protein, but don't add lipids... Just solvate the protein
        # Maybe align along principal axes and then build a cell according to PBC

        # Set PBC starting from diameter and adding distance
        if "cubic".startswith(self.options["pbc"]):
            self.pbcx = self.pbcy = self.pbcz = prot.diam() + self.options["d"]
        elif "rectangular".startswith(self.options["pbc"]):
            self.pbcx, self.pbcy, self.pbcz = vvadd(vvsub(prot.fun(max), prot.fun(min)), self.options["d"])
        else:
            # Rhombic dodecahedron
            self.pbcx = self.pbcy = prot.diam() + self.options["d"]
            self.pbcz = math.sqrt(2) * self.pbcx / 2

        # Possibly override
        self.pbcx = self.pbc_set_x and self.pbc_set_x[0] or self.pbcx
        self.pbcy = self.pbc_set_y and self.pbc_set_y[1] or self.pbcy
        self.pbcz = self.pbc_set_z and self.pbc_set_z[2] or self.pbcz

        # Center coordinates in rectangular brick -- Add solvent next
        if len(self.tm) == 1:
            prot.center((0.5 * self.pbcx, 0.5 * self.pbcy, 0.5 * self.pbcz))

        # Do not set an exclusion range for solvent
        self.options["excl"] = -1
        return prot

    def __process_membrane_no_protein(self):
        """
        Setup membrane without a protein.
        :return: None
        """

        # Set the box -- If there is a disc/hole, add its radius to the distance
        if self.options["disc"]:
            self.pbcx = self.pbcy = self.pbcz = self.options["d"] + 2 * self.options["disc"]
        elif self.options["hole"]:
            self.pbcx = self.pbcy = self.pbcz = self.options["d"] + 2 * self.options["hole"]
        else:
            self.pbcx = self.pbcy = self.pbcz = self.options["d"]

        if "hexagonal".startswith(self.options["pbc"]):
            # Hexagonal prism -- y derived from x directly
            self.pbcy = math.sqrt(3) * self.pbcx / 2
            self.pbcz = self.options["dz"] or self.options["z"] or self.options["d"]
        elif "optimal".startswith(self.options["pbc"]):
            # Rhombic dodecahedron with hexagonal XY plane
            self.pbcy = math.sqrt(3) * self.pbcx / 2
            self.pbcz = math.sqrt(6) * self.options["d"] / 3
        if "rectangular".startswith(self.options["pbc"]):
            self.pbcz = self.options["dz"] or self.options["z"] or self.options["d"]

        # Possibly override
        self.pbcx = self.pbc_set_x and self.pbc_set_x[0] or self.pbcx
        self.pbcy = self.pbc_set_y and self.pbc_set_y[1] or self.pbcy
        self.pbcz = self.pbc_set_z and self.pbc_set_z[2] or self.pbcz

    def write_to_gromacs_or_pdb(self):
        """
        Export membrane in GROMACS or PDB format.
        :return: None
        """
        with open(self.options["o"], "w") as out_file:
            if self.options["o"].endswith(".gro"):
                # Export in GROMACS format
                if self.membrane.atoms:
                    self.title = "INSANE! Membrane UpperLeaflet>" + ":".join(self.lipU) + "=" + ":".join([str(i) for i in self.numU])
                    self.title += " LowerLeaflet>" + ":".join(self.lipL) + "=" + ":".join([str(i) for i in self.numL])

                    if self.protein:
                        self.title = "Protein in " + self.title
                else:
                    self.title = "Insanely solvated protein."

                print(self.title, file=out_file)

                # Print the number of atoms
                print("{:5d}".format(len(self.protein) + len(self.membrane) + len(self.sol)), file=out_file)

                # Print the atoms
                id = 1
                if self.protein:
                    for i in range(len(self.protein)):
                        at, rn, ri = self.protein.atoms[i][:3]
                        x, y, z = self.protein.coord[i]
                        if rn.endswith('.o'):
                            rn = rn[:-2]
                        out_file.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (ri % 1e5, rn, at, id % 1e5, x, y, z))
                        id += 1
                if self.membrane:
                    for i in range(len(self.membrane)):
                        at, rn, ri = self.membrane.atoms[i][:3]
                        x, y, z = self.membrane.coord[i]
                        if rn.endswith('.o'):
                            rn = rn[:-2]
                        out_file.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (ri % 1e5, rn, at, id % 1e5, x, y, z))
                        id += 1
                if self.sol:
                    # Print the solvent
                    print("\n".join([i[0]+"%8.3f%8.3f%8.3f"%i[1] for i in self.sol]), file=out_file)

                # Print the box
                print("%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f" % self.grobox, file=out_file)

            else:
                # Export as PDB format
                if self.membrane.atoms:
                    self.title = "TITLE INSANE! Membrane UpperLeaflet>" + ":".join(self.lipU) + "=" + ":".join([str(i) for i in self.numU])
                    self.title += " LowerLeaflet>" + ":".join(self.lipL) + "=" + ":".join([str(i) for i in self.numL])
                else:
                    self.title = "TITLE Insanely solvated protein."
                print(self.title, file=out_file)

                # Print the box
                print(pdb_box_string(self.box), file=out_file)

                # Print the atoms
                id = 1
                if self.protein:
                    for i in range(len(self.protein)):
                        at, rn, ri = self.protein.atoms[i][:3]
                        x, y, z = self.protein.coord[i]
                        if rn.endswith('.o'):
                            rn = rn[:-2]
                        out_file.write(pdbline % (id % 1e5, at, rn, "", ri % 1e5, '', 10 * x, 10 * y, 10 * z, 0, 0, ''))
                        id += 1
                if self.membrane:
                    for i in range(len(self.membrane)):
                        at, rn, ri = self.membrane.atoms[i][:3]
                        x, y, z = self.membrane.coord[i]
                        if rn.endswith('.o'):
                            rn = rn[:-2]
                        out_file.write(pdbline % (id % 1e5, at, rn, "", ri % 1e5, '', 10 * x, 10 * y, 10 * z, 0, 0, ''))
                        id += 1
                if self.sol:
                    # Print the solvent
                    for i in range(len(self.sol)):
                        ri, rn, at, ai = self.sol[i][0][:5], self.sol[i][0][5:10], self.sol[i][0][10:15], self.sol[i][0][15:20]
                        x, y, z = self.sol[i][1]
                        if rn.endswith('.o'):
                            rn = rn[:-2]
                        out_file.write(pdbline % (
                        id % 1e5, at.strip(), rn.strip(), "", int(ri) % 1e5, '', 10 * x, 10 * y, 10 * z, 0, 0, ''))
                        id += 1

    def __get_topology(self):
        """
        Determine topology.
        :return: Array<tuple>
        """
        topmolecules = []
        for i in self.molecules:
            if i[0].endswith('.o'):
                topmolecules.append(tuple([i[0][:-2]] + list(i[1:])))
            else:
                topmolecules.append(i)
        return topmolecules

    def write_topology(self):
        """
        Write topology to file.
        :return: None
        """
        with open(self.options["p"], "w") as top_file:
            print('#include "martini.itp"\n', file=top_file)
            print('[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%self.title, file=top_file)
            if self.protein:
                print("%-10s %5d"%("Protein", 1), file=top_file)
            print("\n".join("%-10s %7d"%i for i in self.topology), file=top_file)

    @staticmethod
    def _point(y, phi):
        r = math.sqrt(1 - y * y)
        return math.cos(phi) * r, y, math.sin(phi) * r

    @staticmethod
    def points_on_sphere(n):  # Name needs to be refactored
        return [Insane3._point((2. * k + 1) / n - 1, k * 2.3999632297286531) for k in range(n)]


if __name__ == "__main__":
    membrane = Insane3()
