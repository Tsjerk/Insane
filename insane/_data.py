# INSert membrANE
# A simple, versatile tool for building coarse-grained simulation systems
# Copyright (C) 2017  Tsjerk A. Wassenaar and contributors
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
Data tables for building systems.
"""

import numpy as np

__all__ = ['CHARGES', 'SOLVENTS', 'APOLARS']

# Lists for automatic charge determination
CHARGES = {
    "ARG":1, "LYS":1, "HIH": 1, "ASP":-1, "GLU":-1
}

a, b  = np.sqrt(2) / 20, np.sqrt(2) / 60
# Tetrahedral
ct = np.cos(np.pi * 109.47 / 180)
st = np.sin(np.pi * 109.47 / 180)

# The following dict is a data table. Whatever pylint says, it is clearer when
# aligned.
# pylint: disable=bad-whitespace

# Get a set of coordinates for a solvent particle with a given name
# Dictionary of solvents: the keys are the residue names, the values are tuples
# of bead names and coordinates.
# First only the solvents with multiple atoms
SOLVENTS = {
    "PW":       (("W", (-0.07, 0, 0)),  # Polarizable water
                 ("WP", (0.07, 0, 0)),
                 ("WM", (0.07, 0, 0))),
    "BMW":      (("C", (0, 0, 0)),
                 ("Q1", (0.12, 0, 0)),
                 ("Q2", (-0.06, np.cos(np.pi / 6) * 0.12, 0))),  # BMW water
    "SPC":      (("OW", (0, 0, 0)),  # SPC
                 ("HW1", (0.01, 0, 0)),
                 ("HW2", (0.01 * ct, 0.01 * st, 0))),
    "SPCM":     (("OW", (0, 0, 0)),  # Multiscale/Martini SPC
                 ("HW1", (0.01, 0, 0)),
                 ("HW2", (0.01 * ct, 0.01 * st, 0)),
                 ("vW", (0, 0, 0))),
    "FG4W":     (("OW1", (a, a, a)),  # Bundled water
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
    "FG4W-MS":  (("OW1", (a, a, a)),  # Bundled water, multiscaled
                 ("HW11", (a, a - b, a - b)),
                 ("HW12", (a, a + b, a + b)),
                 ("OW2", (a, -a, -a)),
                 ("HW21", (a - b, -a, -a + b)),
                 ("HW22", (a + b, -a, -a -b)),
                 ("OW3", (-a, -a, a)),
                 ("HW31", (-a, -a + b, a - b)),
                 ("HW32", (-a, -a - b, a + b)),
                 ("OW4", (-a, a, -a)),
                 ("HW41", (-a + b, a, -a + b)),
                 ("HW42", (-a - b, a, -a - b)),
                 ("VZ", (0, 0, 0))),
    "GMY":      (("B1", (-0.11, 0,   0)),
                 ("B2", ( 0.05, 0.16, 0)),
                 ("B3", ( 0.05, -0.16, 0))),
    "GLUC":     (("B1", (-0.11, 0,   0)),
                 ("B2", ( 0.05, 0.16, 0)),
                 ("B3", ( 0.05, -0.16, 0))),
    "FRUC":     (("B1", (-0.11, 0,   0)),
                 ("B2", ( 0.05, 0.16, 0)),
                 ("B3", ( 0.05, -0.16, 0))),
    "SUCR":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "MALT":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "CELL":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "KOJI":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "SOPH":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "NIGE":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "LAMI":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
    "TREH":     (("B1", (-0.25, 0.25, 0)),
                 ("B2", (-0.25, 0,   0)),
                 ("B3", (-0.25, -0.25, 0)),
                 ("B4", ( 0.25, 0,   0)),
                 ("B5", ( 0.25, 0.25, 0)),
                 ("B6", ( 0.25, -0.25, 0))),
# Loose aminoacids
    "GLY":      (("BB",  ( 0,    0,   0)), ),
    "ALA":      (("BB",  ( 0,    0,   0)), ),
    "ASN":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "ASP":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "GLU":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "GLN":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "LEU":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "ILE":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "VAL":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "SER":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "THR":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "CYS":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "MET":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "LYS":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "PRO":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "HYP":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", (-0.25, 0,   0))),
    "ARG":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", ( 0,    0,   0)),
                 ("SC2", (-0.25, 0.125, 0))),
    "PHE":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", ( 0,    0,   0)),
                 ("SC2", (-0.25, -0.125, 0)),
                 ("SC3", (-0.25, 0.125, 0))),
    "TYR":      (("BB",  ( 0.25, 0,   0)),
                 ("SC1", ( 0,    0,   0)),
                 ("SC2", (-0.25, -0.125, 0)),
                 ("SC3", (-0.25, 0.125, 0))),
    "TRP":      (("BB",  ( 0.25, 0.125, 0)),
                 ("SC1", ( 0.25, 0,   0)),
                 ("SC2", ( 0,   -0.125, 0)),
                 ("SC3", ( 0,    0.125, 0)),
                 ("SC4", (-0.25, 0,     0))),
    }

# Update the solvents dictionary with single atom ones
for solvent in ["W", "WN", "NA", "CL", "Mg", "K", "BUT"]:
    SOLVENTS[solvent] = ((solvent, (0, 0, 0)), )

# Apolar amino acids nd stuff for orienting proteins in membrane
APOLARS = "ALA CYS PHE ILE LEU MET VAL TRP PLM CLR".split()
