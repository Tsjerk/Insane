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

import sys

import simopt
from simopt import MULTI, MA

from . import core
from .converters import vector, box3d, molspec


# Option list
OPTIONS = simopt.Options([
    #  level opt  attribute      type        num     default    flags    description
        """
    Input/output related options
    """,
        (0, "-f", "solute",      str,         1,        None, MULTI, "Input GRO or PDB file 1: Solute (e.g. Protein)"),
        (0, "-o", "output",      str,         1,        None,    MA, "Output GRO file: Membrane with Protein"),
        (0, "-p", "topology",    str,         1,        None,     0, "Optional rudimentary topology file"),
        """
    Periodic boundary conditions
    If -d is given, set up PBC according to -pbc such that no periodic
    images are closer than the value given. When macromolecules are 
    included -d indicates the distance between any of them, within a 
    unit cell or across PBC. When having only absolute numbers of lipids, 
    the size of the cell will be calculate from the area per lipid.
    """,
        (1, "-pbc", "pbc",       str,         1, "hexagonal",     0, "PBC type: hexagonal, rectangular, square, cubic, optimal or keep"),
        (0, "-d",   "distance",  float,       1,           0,     0, "Distance between periodic images (nm)"),
        (0, "-dz",  "zdistance", float,       1,        None,     0, "Z distance between periodic images (nm)"),
        (2, "-x",   "xvector",   vector,      1,        None,     0, "X dimension or first lattice vector of system (nm)"),
        (2, "-y",   "yvector",   vector,      1,        None,     0, "Y dimension or first lattice vector of system (nm)"),
        (2, "-z",   "zvector",   vector,      1,        None,     0, "Z dimension or first lattice vector of system (nm)"),
        (2, "-box", "box",       box3d,       1,        None,     0, "Box in GRO (3 or 9 floats) or PDB (6 floats) format, comma separated"),
        (0, "-n",   "index",     str,         1,        None,     0, "Index file --- TO BE IMPLEMENTED"),
        """
    Membrane/lipid related options.
    The options -l and -u can be given multiple times. Option -u can be
    used to set the lipid type and abundance for the upper leaflet. Option
    -l sets the type and abundance for the lower leaflet if option -u is
    also given, or for both leaflets if option -u is not given. The
    meaning of the number depends on whether option -d is used to set up
    PBC
    """,
        (0, "-l",    "lower",     molspec,     1,        None, MULTI, "Lipid type and relative abundance (NAME[:#])"),
        (0, "-u",    "upper",     molspec,     1,        None, MULTI, "Lipid type and relative abundance (NAME[:#])"),
        (1, "-a",    "area",      float,       1,        0.60,     0, "Area per lipid (nm*nm)"),
        (1, "-au",   "uparea",    float,       1,        None,     0, "Area per lipid (nm*nm) for upper layer"),
        (1, "-asym", "asymmetry", int,         1,        None,     0, "Membrane asymmetry (number of lipids)"),
        (0, "-hole", "hole",      float,       1,           0,     0, "Make a hole in the membrane with specified radius"),
        (0, "-disc", "disc",      float,       1,        None,     0, "Make a membrane disc with specified radius"),
        (2, "-rand", "randkick",  float,       1,         0.1,     0, "Random kick size (maximum atom displacement)"),
        (2, "-bd",   "beaddist",  float,       1,         0.3,     0, "Bead distance unit for scaling z-coordinates (nm)"),
        """
    Protein related options.
    """,
        (0, "-center", "center",      bool,        0,        None, 0, "Center the protein on z"),
        (9, "-orient", "orient",      bool,        0,        None, 0, "Orient protein in membrane"),
        (1, "-rotate", "rotate",      str,         1,        None, 0, "Rotate protein (random|princ|angle(float)"),
        (9, "-od",     "origriddist", float,       1,         1.0, 0, "Grid spacing for determining orientation"),
        (9, "-op",     "oripower",    float,       1,         4.0, 0, "Hydrophobic ratio power for determining orientation"),
        (2, "-fudge",  "fudge",       float,       1,         0.1, 0, "Fudge factor for allowing lipid-protein overlap"),
        (1, "-ring",   "inside",      bool,        0,        None, 0, "Put lipids inside the protein"),
        (1, "-dm",     "memshift",    float,       1,           0, 0, "Shift protein with respect to membrane"),
        """
    Solvent related options.
    """,
        (0, "-sol",    "solvent",     molspec,     1,        None, MULTI, "Solvent type and relative abundance (NAME[:#])"),
        (1, "-sold",   "soldiam",     float,       1,         0.5,     0, "Solvent diameter"),
        (1, "-solr",   "solrandom",   float,       1,         0.1,     0, "Solvent random kick"),
        (2, "-excl",   "solexcl",     float,       1,         1.5,     0, "Exclusion range (nm) for solvent addition relative to membrane center"),
        """
    Salt related options.
    """,
        (0, "-salt",   "salt",        str,         1,        None,     0, "Salt concentration"),
        (1, "-charge", "charge",      str,         1,      "auto",     0, "Charge of system. Set to auto to infer from residue names"),
        """
    Define additional lipid types (same format as in lipid-martini-itp-v01.py)
    """,
        (1, "-alname",   "lipnames",  str,  1,  None, MULTI, "Additional lipid name, x4 letter"),
        (1, "-alhead",   "lipheads",  str,  1,  None, MULTI, "Additional lipid head specification string"),
        (1, "-allink",   "liplinks",  str,  1,  None, MULTI, "Additional lipid linker specification string"),
        (1, "-altail",   "liptails",  str,  1,  None, MULTI, "Additional lipid tail specification string"),
        (1, "-alcharge", "lipcharge", str,  1,  None, MULTI, "Additional lipid charge"),
        (0, "-m",        "molfile",   str,  1,  None, MULTI, "Read molecule definitions from file"),
        ])


def main(argv):
    ## TEMPORARY ---
    # Exception to be defined in insane
    class InsaneBuildException(BaseException): pass
    ## <---

    ## OPTIONS
    # Parse options
    try:
        options = OPTIONS.parse(argv[1:])
    except simopt.SimoptHelp:
        print(OPTIONS.help(argv[1:]))
        return 0
    except simopt.MissingMandatoryError as e:
        print(e)
        return 3
    except simopt.Usage as e:
        print(e)
        return 1

    ## WORK
    try:
        system = core.insane(**options)
    except InsaneBuildException as e:
        print(e)
        return 2

    ## OUTPUT
    # Build atom list
    # Build topology
    # Build index



    (molecules,
     protein,
     membrane,
     solvent,
     lipids,
     box) = core.old_main(**options)

    title = core.system_title(membrane, protein, lipids)
    atoms = protein + membrane + solvent

    core.write_summary(protein, membrane, solvent)
    core.write_structure(
        output=options['output'],
        title=title,
        atoms=atoms,
        box=box,
    )
    core.write_top(options['topology'], molecules, title)

    return 0


def cli():
    sys.exit(main(sys.argv))
