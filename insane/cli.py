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
from simopt import MULTI

from . import core
from .converters import vector, box3d, molspec


# Option list
OPTIONS = simopt.Options([
    #   opt          attribute        type  num     default    flags    description
    #   option           type number default description
        """
    Input/output related options
    """,
        ("-f", "solute",      str,         1,        None, MULTI, "Input GRO or PDB file 1: Solute (e.g. Protein)"),
        ("-o", "output",      str,         1,        None,     0, "Output GRO file: Membrane with Protein"),
        ("-p", "topology",    str,         1,        None,     0, "Optional rudimentary topology file"),
        """
    Periodic boundary conditions
    If -d is given, set up PBC according to -pbc such that no periodic
    images are closer than the value given.  This will make the numbers
    provided for lipids be interpreted as relative numbers. If -d is
    omitted, those numbers are interpreted as absolute numbers, and the
    PBC are set to fit the given number of lipids in.
    """,
        ("-pbc", "pbc",       str,         1, "hexagonal",     0, "PBC type: hexagonal, rectangular, square, cubic, optimal or keep"),
        ("-d",   "distance",  float,       1,           0,     0, "Distance between periodic images (nm)"),
        ("-dz",  "zdistance", float,       1,        None,     0, "Z distance between periodic images (nm)"),
        ("-x",   "xvector",   vector,      1,        None,     0, "X dimension or first lattice vector of system (nm)"),
        ("-y",   "yvector",   vector,      1,        None,     0, "Y dimension or first lattice vector of system (nm)"),
        ("-z",   "zvector",   vector,      1,        None,     0, "Z dimension or first lattice vector of system (nm)"),
        ("-box", "box",       box3d,       1,        None,     0, "Box in GRO (3 or 9 floats) or PDB (6 floats) format, comma separated"),
        ("-n",   "index",     str,         1,        None,     0, "Index file --- TO BE IMPLEMENTED"),
        """
    Membrane/lipid related options.
    The options -l and -u can be given multiple times. Option -u can be
    used to set the lipid type and abundance for the upper leaflet. Option
    -l sets the type and abundance for the lower leaflet if option -u is
    also given, or for both leaflets if option -u is not given. The
    meaning of the number depends on whether option -d is used to set up
    PBC
    """,
        ("-l",    "lower",     molspec,     1,        None, MULTI, "Lipid type and relative abundance (NAME[:#])"),
        ("-u",    "upper",     molspec,     1,        None, MULTI, "Lipid type and relative abundance (NAME[:#])"),
        ("-a",    "area",      float,       1,        0.60,     0, "Area per lipid (nm*nm)"),
        ("-au",   "uparea",    float,       1,        None,     0, "Area per lipid (nm*nm) for upper layer"),
        ("-asym", "asymmetry", int,         1,        None,     0, "Membrane asymmetry (number of lipids)"),
        ("-hole", "hole",      float,       1,           0,     0, "Make a hole in the membrane with specified radius"),
        ("-disc", "disc",      float,       1,        None,     0, "Make a membrane disc with specified radius"),
        ("-rand", "randkick",  float,       1,         0.1,     0, "Random kick size (maximum atom displacement)"),
        ("-bd",   "beaddist",  float,       1,         0.3,     0, "Bead distance unit for scaling z-coordinates (nm)"),
        """
    Protein related options.
    """,
        ("-center", "center",      bool,        0,        None, 0, "Center the protein on z"),
        ("-orient", "orient",      bool,        0,        None, 0, "Orient protein in membrane"),
        ("-rotate", "rotate",      str,         1,        None, 0, "Rotate protein (random|princ|angle(float)"),
        ("-od",     "origriddist", float,       1,         1.0, 0, "Grid spacing for determining orientation"),
        ("-op",     "oripower",    float,       1,         4.0, 0, "Hydrophobic ratio power for determining orientation"),
        ("-fudge",  "fudge",       float,       1,         0.1, 0, "Fudge factor for allowing lipid-protein overlap"),
        ("-ring",   "inside",      bool,        0,        None, 0, "Put lipids inside the protein"),
        ("-dm",     "memshift",    float,       1,           0, 0, "Shift protein with respect to membrane"),
        """
    Solvent related options.
    """,
        ("-sol",    "solvent",     molspec,     1,        None, MULTI, "Solvent type and relative abundance (NAME[:#])"),
        ("-sold",   "soldiam",     float,       1,         0.5,     0, "Solvent diameter"),
        ("-solr",   "solrandom",   float,       1,         0.1,     0, "Solvent random kick"),
        ("-excl",   "solexcl",     float,       1,         1.5,     0, "Exclusion range (nm) for solvent addition relative to membrane center"),
        """
    Salt related options.
    """,
        ("-salt",   "salt",        str,         1,        None,     0, "Salt concentration"),
        ("-charge", "charge",      str,         1,      "auto",     0, "Charge of system. Set to auto to infer from residue names"),
        """
    Define additional lipid types (same format as in lipid-martini-itp-v01.py)
    """,
        ("-alname",   "lipnames",  str,  1,  None, MULTI, "Additional lipid name, x4 letter"),
        ("-alhead",   "lipheads",  str,  1,  None, MULTI, "Additional lipid head specification string"),
        ("-allink",   "liplinks",  str,  1,  None, MULTI, "Additional lipid linker specification string"),
        ("-altail",   "liptails",  str,  1,  None, MULTI, "Additional lipid tail specification string"),
        ("-alcharge", "lipcharge", str,  1,  None, MULTI, "Additional lipid charge"),
        ("-m",        "molfile",   str,  1,  None, MULTI, "Read molecule definitions from file"),
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
     lipU, lipL,
     numU, numL,
     box) = core.old_main(argv, options)
    core.write_all(
        output=options['output'],
        topology=options['topology'],
        molecules=molecules,
        protein=protein,
        membrane=membrane,
        solvent=solvent,
        lipU=lipU,
        lipL=lipL,
        numU=numU,
        numL=numL,
        box=box,
    )

    return 0


def cli():
    sys.exit(main(sys.argv))
