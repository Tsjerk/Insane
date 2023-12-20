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

from collections.abc import MutableMapping
import math
import os

from . import utils

__all__ = ['Lipid', 'Lipid_List', 'get_lipids']


# Lipid data file
LIPID_FILE = 'lipids.dat'

# Define supported lipid head beads. Name mapped to atom name
HEADBEADS = {
    "C":  "NC3", # NC3 = Choline
    "E":  "NH3", # NH3 = Ethanolamine
    "G":  "GL0", # GL0 = Glycerol
    "S":  "CNO", # CNO = Serine
    "P":  "PO4", # PO4 = Phosphate
    "PS1":"PS1", # PS1 is bead one of two bead PS represents a COO group
    "PS2":"PS2", # PS2 is bead two of two bead PS represents a NH3 group
    "COH":"COH", # COH = Cappding bead for top of diacylglycerols and ceramides
}

# Define supported lipid link beads. Name mapped to atom name
LINKBEADS = {
    "G":  "GL",  # Glycerol
    "A":  "AM",  # Amide (Ceramide/Sphingomyelin)
    "O":  "OH",  # Amide (Ceramide/Sphingomyelin), in Martini 3 old AM1 is called OH1 
}


class Lipid:
    """Lipid structure"""

    def __init__(self, **kwargs):
        self.name      = kwargs.get("name")
        self.head      = kwargs.get("head")
        self.link      = kwargs.get("link")
        self.tail      = kwargs.get("tail")
        self.beads     = kwargs.get("beads")
        if type(self.beads) == str:
            self.beads = self.beads.split()
        self.charge    = kwargs.get("charge")
        self.template  = kwargs.get("template")
        self.area      = kwargs.get("area")
        self.diam      = kwargs.get("diam", math.sqrt(kwargs.get("area", 0)))
        self.coords    = None
        if kwargs.get("string"):
            self.parse(kwargs["string"])

    def parse(self, string):
        """
        Parse lipid definition from string:

            alhead=C P, allink=A A, altail=TCC CCCC, alname=DPSM, charge=0.0
        """
        fields = [i.split("=") for i in string.split(', ')]
        for what, val in fields:
            what = what.strip()
            val  = val.split()
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
        # Not been impleneted yet 
        #if self.charge is None:
        #    # Infer charge from head groups
        #    self.charge = sum([headgroup_charges[bead] for bead in self.head])

    def build(self, **kwargs):
        """Build/return a list of [(bead, x, y, z), ...]"""

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
                    beads = [HEADBEADS[i] for i in self.head]
                    beads.extend([LINKBEADS[n] + str(i + 1)
                                  for i, n in enumerate(self.link)])
                    for i, t in enumerate(self.tail):
                        beads.extend([n + chr(65 + i) + str(j + 1)
                                      for j, n in enumerate(t)])

                taillength = max([0]+[len(i) for i in self.tail])
                length = len(self.head)+taillength

                # Add the pseudocoordinates for the head
                rl     = range(len(self.head))
                struc  = [(0, 0, length-i) for i in rl]

                # Add the linkers
                rl     = range(len(self.link))
                struc.extend([(i%2, i//2, taillength) for i in rl ])

                # Add the tails
                for j, tail in enumerate(self.tail):
                    rl = range(len(tail))
                    struc.extend([(j%2, j//2, taillength-1-i) for i in rl])

                mx, my, mz = [ (max(i)+min(i))/2 for i in zip(*struc) ]
                self.coords = [
                    [i, 0.25 * (x - mx), 0.25 * (y - my), z]
                    for i, (x, y, z) in zip(beads, struc)
                ]

        # Scale the x/y based on the lipid's APL - diameter is less than sqrt(APL)
        diam   = kwargs.get("diam", self.diam)
        radius = diam*0.45
        minmax = [ (min(i), max(i)) for i in list(zip(*self.coords))[1:] ]
        mx, my, mz = [ sum(i)/2. for i in minmax ]
        scale  = radius/math.sqrt((minmax[0][0]-mx)**2 + (minmax[1][0]-my)**2)

        for i in self.coords:
            i[1] = scale*(i[1]-mx)
            i[2] = scale*(i[2]-my)
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


class Lipid_List(MutableMapping):
    """Container class for lipid definitions"""

    def __init__(self):
        self.store = dict()
        self.last  = None

    def __getitem__(self, key):
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

    def add_from_stream(self, stream):
        current_lipid = None
        for line in stream:
            if line[0] == ';':
                line = line[1:]
            if line.startswith("@INSANE"):
                current_lipid = Lipid(string=line[len('@INSANE'):].strip())
                self[current_lipid.name] = current_lipid
            elif line.startswith("@BEADS"):
                beads = line.split()[1:]
                current_lipid.beads = beads

    def add_from_file(self, path):
        with open(path) as infile:
            self.add_from_stream(infile)

    def add_from_files(self, multi_path):
        for path in multi_path:
            self.add_from_file(path)

    def add_from_def(self, usrnames, usrheads, usrlinks, usrtails, usrcharges):
        if not usrcharges:
            usrcharges = [None] * len(usrnames)
        lipids = zip(usrnames, usrheads, usrlinks, usrtails, usrcharges)
        for name, head, link, tail, charge in lipids:
            heads = head.replace(".", " ").split()
            linkers = link.replace(".", " ").split()
            tails = tail.replace(".", " ").split()
            charge = None if charge is None else float(charge)
            self[name] = Lipid(name=name,
                               head=heads,
                               link=linkers,
                               tail=tails,
                               charge=charge)


def read_lipids(lipfile):
    lipids = Lipid_List()
    x, y, z = None, None, None
    for line in lipfile:
        stripped = line.strip()
        if (not stripped) or stripped.startswith((';','#')):
            # Comment or empty line
            continue
        elif stripped.startswith('['):
            # Moleculetype tag
            x, y, z = None, None, None
            # moltype = stripped[2:stripped.find(']')].strip()
        elif x is None:
            x = [float(i) for i in line.split() ]
        elif y is None:
            y = [float(i) for i in line.split() ]
        elif z is None:
            z = [float(i) for i in line.split() ]
        else:
            splitted = line.split()
            name = splitted.pop(0)
            charge = splitted.pop(0)
            lipids[name] = Lipid(
                name=name, 
                charge=charge, 
                beads=splitted, 
                template=zip(x, y, z)
            )
    return lipids


def get_lipids():
    lipid_stream = utils.iter_resource(LIPID_FILE)
    lipids = read_lipids(lipid_stream)
    return lipids
