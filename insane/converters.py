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

import math

from .constants import d2r


# These functions are not to stay here... -TAW
def vector(v):
    if type(v) == str and "," in v:
        return [float(i) for i in v.split(",")]
    return float(v)


def pdbBoxRead(a):
    # Convert a PDB CRYST1 entry to a lattice definition.
    # Convert from Angstrom to nanometer
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac) , math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [0.1*fa, 0, 0, 0.1*fb*cg, 0.1*fb*sg, 0, wx, wy, wz]


def box3d(a):
    x = [ float(i) for i in a.split(",") ] + 6*[0]
    if len(x) == 12: # PDB format
        return pdbBoxRead("CRYST1 "+" ".join([str(i) for i in x]))
    elif len(x) == 7:
        return x[0], x[0], x[0], 0, 0, 0, 0, 0, 0
    else:            # GRO format
        return x[0], x[3], x[4], x[5], x[1], x[6], x[7], x[8], x[2]


def molspec(x):
    """
    Parse a string for a lipid or a solvent as given on the command line
    (MOLECULE[=NUMBER|:NUMBER]); where `=NUMBER` sets an absolute number of the
    molecule, and `:NUMBER` sets a relative number of it.
    If both absolute and relative number are set False, then relative count is 1.
    """
    lip = x.split(":")
    abn = lip[0].split("=")
    names = abn[0]
    if len(abn) > 1:
        nrel = 0
        nabs = int(abn[1])
    else:
        nabs = 0
        if len(lip) > 1:
            nrel = float(lip[1])
        else:
            nrel = 1
    return abn[0], nabs, nrel
