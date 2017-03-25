#!/usr/bin/env python

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
INSANE: A versatile tool for building membranes and/or solvent with proteins.

... Someone ought to write a more extensive docstring here...
"""


__authors__ = ["Tsjerk A. Wassenaar", "Jonathan Barnoud"]
__version__ = "1.0-dev"
__year__    = "2017"

version   = "20160625.15.TAW"
previous  = "20160104.18.TAW"


import os
import sys
import math
import random
import collections

from . import linalg
from .converters import *
from .constants import d2r
from .data import (lipidsa, lipidsx, lipidsy, lipidsz, headbeads, linkbeads,
                   solventParticles, charges, apolar)


# Set the random seed.
# The seed is set to an arbitary value set in the INSANE_SEED environment
# variable. If the environment variable is not set, then the system time is
# used to set the seed.
random.seed(os.environ.get('INSANE_SEED', None))


## PRIVATE PARTS FROM THIS POINT ON ##

def mean(a):
    return sum(a)/len(a)

def isPDBAtom(l):
    return l.startswith("ATOM") or l.startswith("HETATM")

def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>   atom name,   res name,     res id, chain,       x,            y,             z
    return (a[12:16], a[17:20], int(a[22:26]), a[21], float(a[30:38])/10, float(a[38:46])/10, float(a[46:54])/10)


# Reformatting of lines in structure file
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"
pdbline = "ATOM  %5i  %-3s %4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s  \n"


def pdbBoxString(box):
    # Box vectors
    u, v, w  = box

    # Box vector lengths
    nu, nv, nw = [math.sqrt(linalg.norm2(i)) for i in (u, v, w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(linalg.cos_angle(v, w))/d2r
    beta  = nu*nw == 0 and 90 or math.acos(linalg.cos_angle(u, w))/d2r
    gamma = nu*nv == 0 and 90 or math.acos(linalg.cos_angle(u, v))/d2r

    return pdbBoxLine % (10*linalg.norm(u),
                         10*linalg.norm(v),
                         10*linalg.norm(w),
                         alpha, beta, gamma)


def groAtom(a):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>   atom name,   res name,     res id, chain,       x,          y,          z
    return (a[10:15], a[5:10],   int(a[:5]), " ", float(a[20:28]), float(a[28:36]), float(a[36:44]))

def groBoxRead(a):
    b = [float(i) for i in a.split()] + 6*[0] # Padding for rectangular boxes
    return b[0], b[3], b[4], b[5], b[1], b[6], b[7], b[8], b[2]


class Structure:
    def __init__(self, filename=None):
        self.title   = ""
        self.atoms   = []
        self.coord   = []
        self.rest    = []
        self.box     = []
        self._center = None

        if filename:
            lines = open(filename).readlines()
            # Try extracting PDB atom/hetatm definitions
            self.rest   = []
            self.atoms  = [pdbAtom(i) for i in lines if isPDBAtom(i) or self.rest.append(i)]
            if self.atoms:
                # This must be a PDB file
                self.title = "THIS IS INSANE!\n"
                for i in self.rest:
                    if i.startswith("TITLE"):
                        self.title = i
                self.box   = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                for i in self.rest:
                    if i.startswith("CRYST1"):
                        self.box = pdbBoxRead(i)
            else:
                # This should be a GRO file
                self.atoms = [groAtom(i) for i in lines[2:-1]]
                self.rest  = [lines[0], lines[1], lines[-1]]
                self.box   = groBoxRead(lines[-1])
                self.title = lines[0]
            self.coord = [i[4:7] for i in self.atoms]
            self.center()

    def __nonzero__(self):
        return bool(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __iadd__(self, s):
        for i in range(len(self)):
            self.coord[i] = linalg.vvadd(self.coord[i], s)
        return self

    def center(self, other=None):
        if not self._center:
            self._center = [ sum(i)/len(i) for i in zip(*self.coord)]
        if other:
            s = linalg.vvsub(other, self._center)
            for i in range(len(self)):
                self.coord[i] = linalg.vvadd(self.coord[i], s)
            self._center = other
            return s # return the shift
        return self._center

    def diam(self):
        if self._center != (0, 0, 0):
            self.center((0, 0, 0))
        return 2*math.sqrt(max([i*i+j*j+k*k for i, j, k in self.coord]))

    def diamxy(self):
        if self._center != (0, 0, 0):
            self.center((0, 0, 0))
        return 2*math.sqrt(max([i*i+j*j for i, j, k in self.coord]))

    def fun(self, fn):
        return [fn(i) for i in zip(*self.coord)]

    def orient(self, d, pw):
        # Determine grid size
        mx, my, mz = self.fun(min)
        rx, ry, rz = self.fun(lambda x: max(x)-min(x)+1e-8)

        # Number of grid cells
        nx, ny, nz = int(rx/d+0.5), int(ry/d+0.5), int(rz/d+0.5)

        # Initialize grids
        atom     = [[[0 for i in range(nz+2)]
                     for j in range(ny+2)] for k in range(nx+2)]
        phobic   = [[[0 for i in range(nz+2)]
                     for j in range(ny+2)] for k in range(nx+2)]
        surface  = []
        for i, (ix, iy, iz) in zip(self.atoms, self.coord):
            if i[1] != "DUM":
                jx, jy, jz = (int(nx*(ix-mx)/rx),
                              int(ny*(iy-my)/ry),
                              int(nz*(iz-mz)/rz))
                atom[jx][jy][jz]   += 1
                phobic[jx][jy][jz] += (i[1].strip() in apolar)

        # Determine average density
        occupd = sum([bool(k) for i in atom for j in i for k in j])
        avdens = float(sum([sum(j) for i in atom for j in i]))/occupd

        #cgofile  = open('density.cgo', "w")
        #cgofile.write('[\n')
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if atom[i][j][k] > 0.1*avdens:
                        # Check the neighbouring cells;
                        # if one of them is not occupied, count cell as surface
                        if not (atom[i-1][j][k] and atom[i+1][j][k] and
                                atom[i][j-1][k] and atom[i][j+1][k] and
                                atom[i][j][k-1] and atom[i][j][k+1]):
                            sx, sy, sz = (mx+rx*(i+0.5)/nx,
                                          my+ry*(j+0.5)/ny,
                                          mz+rz*(k+0.5)/nz)
                            sw       = (2.0*phobic[i][j][k]/atom[i][j][k])**pw
                            surface.append((sx, sy, sz, sw))
                            #cgofile.write("    7.0, %f, %f, %f, %f, \n"%(10*sx, 10*sy, 10*sz, 0.25*sw))
        #cgofile.write(']\n')
        #cgofile.close()

        sx, sy, sz, w = zip(*surface)
        W             = 1.0/sum(w)

        # Weighted center of apolar region; has to go to (0, 0, 0)
        sxm, sym, szm   = [sum(p)*W
                         for p in zip(*[(m*i, m*j, m*k)
                                        for m, i, j, k in zip(w, sx, sy, sz)])]

        # Place apolar center at origin
        self.center((-sxm, -sym, -szm))
        sx, sy, sz    = zip(*[(i-sxm, j-sym, k-szm) for i, j, k in zip(sx, sy, sz)])

        # Determine weighted deviations from centers
        dx, dy, dz      = zip(*[(m*i, m*j, m*k) for m, i, j, k in zip(w, sx, sy, sz)])

        # Covariance matrix for surface
        xx, yy, zz, xy, yz, zx = [sum(p)*W
                             for p in zip(*[(i*i, j*j, k*k, i*j, j*k, k*i)
                                            for i, j, k in zip(dx, dy, dz)])]

        # PCA: u, v, w are a rotation matrix
        (ux, uy, uz), (vx, vy, vz), (wx, wy, wz), r = linalg.mijn_eigen_sym_3x3(xx, yy, zz, xy, zx, yz)

        # Rotate the coordinates
        self.coord = [(ux*i+uy*j+uz*k, vx*i+vy*j+vz*k, wx*i+wy*j+wz*k)
                      for i, j, k in self.coord]

    def rotate_princ(self):
        x, y, z = zip(*self.coord)

        # The rotation matrix in the plane equals the transpose
        # of the matrix of eigenvectors from the 2x2 covariance
        # matrix of the positions.
        # For numerical stability we do
        # d_i     = x_i - x_0
        # mean(x) = x_0 + sum(d_i)/N =
        # var(x)  = sum((d_i - mean(d))**2)/(N-1)
        xy        = ssd(x, y)
        if xy != 0:
            xx     = ssd(x, x)
            yy     = ssd(y, y)

            # The eigenvalues are the roots of the 2nd order
            # characteristic polynomial, with the coefficients
            # equal to the trace and the determinant of the
            # matrix.
            t,  d  = xx+yy, xx*yy - xy*xy
            # The two eigenvectors form a 2D rotation matrix
            # R = ((cos, sin), (-sin, cos)), which means that
            # the second eigenvector follows directly from
            # the first. We thus only need to determine one.
            l1     = t/2 + math.sqrt(0.25*t*t-d)

            ux, uy = l1-yy, xy
            lu     = math.sqrt(ux*ux+uy*uy)

            ux    /=  lu
            uy    /=  lu

            # Finally we rotate the system in the plane by
            # matrix multiplication with the transpose of
            # the matrix of eigenvectors
            self.coord = [(ux*i+uy*j, ux*j-uy*i, k) for i, j, k in zip(x, y, z)]

    def rotate_random(self):
        ux   = math.cos(random.random()*2*math.pi)
        uy   = math.sqrt(1-ux*ux)
        self.coord = [(ux*i+uy*j, ux*j-uy*i, k) for i, j, k in self.coord]

    def rotate_degrees(self, angle):
        ux   = math.cos(angle*math.pi/180.)
        uy   = math.sin(angle*math.pi/180.)
        self.coord = [(ux*i+uy*j, ux*j-uy*i, k) for i, j, k in self.coord]


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
        if self.charge is None:
            # Infer charge from head groups
            self.charge = sum([headgroup_charges[bead] for bead in self.head])

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
                    beads = [ headbeads[i] for i in self.head ]
                    beads.extend([ linkbeads[n]+str(i+1) for i, n in enumerate(self.link)])
                    for i, t in enumerate(self.tail):
                        beads.extend([n+chr(65+i)+str(j+1) for j, n in enumerate(t)])

                taillength = max([0]+[len(i) for i in self.tail])
                length = len(self.head)+taillength

                # Add the pseudocoordinates for the head
                rl     = range(len(self.head))
                struc  = [(0, 0, length-i) for i in rl]

                # Add the linkers
                rl     = range(len(self.link))
                struc.extend([(i%2, i/2, taillength) for i in rl ])

                # Add the tails
                for j, tail in enumerate(self.tail):
                    rl = range(len(tail))
                    struc.extend([(j%2, j/2, taillength-1-i) for i in rl])

                mx, my, mz = [ (max(i)+min(i))/2 for i in zip(*struc) ]
                self.coords = [[i, 0.25*(x-mx), 0.25*(y-my), z] for i, (x, y, z) in zip(beads, struc)]

        # Scale the x/y based on the lipid's APL - diameter is less than sqrt(APL)
        diam   = kwargs.get("diam", self.diam)
        radius = diam*0.45
        minmax = [ (min(i), max(i)) for i in zip(*self.coords)[1:] ]
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


class Lipid_List(collections.MutableMapping):
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


# Mean of deviations from initial value
def meand(v):
    return sum([i-v[0] for i in v])/len(v)

# Sum of squares/crossproducts of deviations
def ssd(u, v):
    return sum([(i-u[0])*(j-v[0]) for i, j in zip(u, v)])/(len(u)-1)

# Parse a string for a lipid as given on the command line (LIPID[:NUMBER])
def parse_mol(x):
    l = x.split(":")
    return l[0], len(l) == 1 and 1 or float(l[1])


# Very simple option class
class Option:
    def __init__(self, func=str, num=1, default=None, description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self):
        return self.value != None
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self, v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]


def old_main(argv, options):

    lipL      = options["lower"]
    lipU      = options["upper"]
    solv      = options["solvent"]
    usrlip    = Lipid_List()
    mollist   = options["molfile"]
    usrnames  = options["lipnames"]
    usrheads  = options["lipheads"]
    usrlinks  = options["liplinks"]
    usrtails  = options["liptails"]
    usrcharg  = options["lipcharge"]

    numU = []
    numL = []

    # Description
    desc = ""

    # Read in the structures (if any)
    tm    = [ Structure(i) for i in options["solute"] ]


    ## I. STRUCTURES

    ## a. LIPIDS

    liplist = Lipid_List()


    ## ==> LIPID  BOOKKEEPING:
    # import lipids
    # lipid_list = lipids.get_list()
    # lipid_list.add_from_file(*options["mollist"])
    # lipid_list.add_from_def(*zip(options["lipnames"], options["lipheads"], options["liplinks"], options["liptails"], options["lipcharge']))


    # First add internal lipids
    for name, lip in lipidsa.items():
        moltype  = lip[0]
        template = zip(lipidsx[moltype], lipidsy[moltype], lipidsz[moltype])
        liplist[name] = Lipid(name=name, beads=lip[1], template=template)

    # Then add lipids from file
    for filename in options["molfile"]:
        stuff = open(filename).read().split("@INSANE")
        for group in stuff[1:]:
            lines  = group.split("\n")
            lipdef = lines.pop(0)
            beads  = None
            for line in lines:
                if line.startswith('[') or not line.strip():
                    break
                if "@BEADS" in line:
                    beads = line.split("@BEADS")[1].split()
            lip = Lipid(string=lipdef, beads=beads)
            liplist[lip.name] = lip

    # Last, add lipids from command line
    for name, head, link, tail in zip(usrnames, usrheads, usrlinks, usrtails):
        heads   = head.replace(".", " ").split()
        linkers = link.replace(".", " ").split()
        tails   = tail.replace(".", " ").split()
        liplist[name] = Lipid(name=name, head=heads, link=linkers, tail=tails)

    # <=== END OF LIPID BOOKKEEPING


    # Periodic boundary conditions

    # option -box overrides everything
    if options["box"]:
        options["xvector"] = options["box"][:3]
        options["yvector"] = options["box"][3:6]
        options["zvector"] = options["box"][6:]

    # option -pbc keep really overrides everything
    if options["pbc"] == "keep" and tm:
        options["xvector"] = tm[0].box[:3]
        options["yvector"] = tm[0].box[3:6]
        options["zvector"] = tm[0].box[6:]

    # options -x, -y, -z take precedence over automatic determination
    pbcSetX = 0
    if type(options["xvector"]) in (list, tuple):
        pbcSetX = options["xvector"]
    elif options["xvector"]:
        pbcSetX = [options["xvector"], 0, 0]

    pbcSetY = 0
    if type(options["yvector"]) in (list, tuple):
        pbcSetY = options["yvector"]
    elif options["yvector"]:
        pbcSetY = [0, options["yvector"], 0]

    pbcSetZ = 0
    if type(options["zvector"]) in (list, tuple):
        pbcSetZ = options["zvector"]
    elif options["zvector"]:
        pbcSetZ = [0, 0, options["zvector"]]


    lo_lipd  = math.sqrt(options["area"])
    if options["uparea"]:
        up_lipd = math.sqrt(options["uparea"])
    else:
        up_lipd = lo_lipd


    ################
    ## I. PROTEIN ##
    ################


    protein  = Structure()
    prot     = []
    xshifts  = [0] # Shift in x direction per protein

    ## A. NO PROTEIN ---
    if not tm:

        resi = 0

        # Set the box -- If there is a disc/hole, add its radius to the distance
        if options["disc"]:
            pbcx = pbcy = pbcz = options["distance"] + 2*options["disc"]
        elif options["hole"]:
            pbcx = pbcy = pbcz = options["distance"] + 2*options["hole"]
        else:
            pbcx = pbcy = pbcz = options["distance"]

        if "hexagonal".startswith(options["pbc"]):
            # Hexagonal prism -- y derived from x directly
            pbcy = math.sqrt(3)*pbcx/2
            pbcz = (options["zdistance"]
                    or options["zvector"]
                    or options["distance"])
        elif "optimal".startswith(options["pbc"]):
            # Rhombic dodecahedron with hexagonal XY plane
            pbcy = math.sqrt(3)*pbcx/2
            pbcz = math.sqrt(6)*options["distance"]/3
        if "rectangular".startswith(options["pbc"]):
            pbcz = (options["zdistance"]
                    or options["zvector"]
                    or options["distance"])

        # Possibly override
        pbcx = pbcSetX and pbcSetX[0] or pbcx
        pbcy = pbcSetY and pbcSetY[1] or pbcy
        pbcz = pbcSetZ and pbcSetZ[2] or pbcz


    ## B. PROTEIN ---
    else:

        for prot in tm:


            ## a. NO MEMBRANE --
            if not lipL:

                # A protein, but don't add lipids... Just solvate the protein
                # Maybe align along principal axes and then build a cell according to PBC

                # Set PBC starting from diameter and adding distance
                if "cubic".startswith(options["pbc"]):
                    pbcx = pbcy = pbcz = prot.diam() + options["distance"]
                elif "rectangular".startswith(options["pbc"]):
                    pbcx, pbcy, pbcz = linalg.vvadd(linalg.vvsub(prot.fun(max), prot.fun(min)), options["distance"])
                else:
                    # Rhombic dodecahedron
                    pbcx = pbcy = prot.diam()+options["distance"]
                    pbcz = math.sqrt(2)*pbcx/2

                # Possibly override
                pbcx = pbcSetX and pbcSetX[0] or pbcx
                pbcy = pbcSetY and pbcSetY[1] or pbcy
                pbcz = pbcSetZ and pbcSetZ[2] or pbcz

                # Center coordinates in rectangular brick -- Add solvent next
                if len(tm) == 1:
                    prot.center((0.5*pbcx, 0.5*pbcy, 0.5*pbcz))

                # Do not set an exclusion range for solvent
                options["solexcl"] = -1


            ## b. PROTEIN AND MEMBRANE --
            else:

                # Have to build a membrane around the protein.
                # So first put the protein in properly.

                # Center the protein and store the shift
                shift = prot.center((0, 0, 0))

                ## 1. Orient with respect to membrane
                # Orient the protein according to the TM region, if requested
                # This doesn't actually work very well...
                if options["orient"]:
                    # Grid spacing (nm)
                    d  = options["origriddist"]
                    pw = options["oripower"]

                    prot.orient(d, pw)

                ## 4. Orient the protein in the xy-plane
                ## i. According to principal axes and unit cell
                if options["rotate"] == "princ":
                    prot.rotate_princ()
                ## ii. Randomly
                elif options["rotate"] == "random":
                    prot.rotate_random()

                ## iii. Specifically
                elif options["rotate"]:
                    prot.rotate_degrees(float(options["rotate"]))

                ## 5. Determine the minimum and maximum x and y of the protein
                pmin, pmax = prot.fun(min), prot.fun(max)
                prng       = (pmax[0]-pmin[0], pmax[1]-pmin[1], pmax[2]-pmin[2])
                center     = (0.5*(pmin[0]+pmax[0]), 0.5*(pmin[1]+pmax[1]))


                # Set the z-dimension
                pbcz  = pbcSetZ and pbcSetZ[2]
                # If it is not set, set pbcz to the dimension of the protein
                pbcz  = pbcz or prng[2]
                pbcz += options["zdistance"] or options["distance"] or 0


                # At this point we should shift the subsequent proteins such
                # that they end up at the specified distance, in case we have
                # a number of them to do
                # y-shift is always -ycenter
                # x-shift is -xmin+distance+xmax(current)
                xshft, yshft = (xshifts[-1]-pmin[0]+(options["distance"] or 0),
                                -center[1])
                xshifts.append(xshifts[-1]+pmax[0]+(options["distance"] or 0))


                ## 6. Set box (brick) dimensions
                if options["disc"]:
                    pbcx = options["distance"] + 2*options["disc"]
                    if ("square".startswith(options["pbc"]) or
                        "rectangular".startswith(options["pbc"])):
                        pbcy = pbcx
                    else:
                        pbcy  = math.cos(math.pi/6)*pbcx
                else:
                    pbcx = (options["distance"] or 0) + prng[0]
                    if "square".startswith(options["pbc"]):
                        pbcy = pbcx
                    elif "rectangular".startswith(options["pbc"]):
                        pbcy = options["distance"] + prng[1]
                    else:
                        # This goes for a hexagonal cell as well as for the optimal arrangement
                        # The latter is hexagonal in the membrane plane anyway...
                        pbcy  = math.cos(math.pi/6)*pbcx


                ## 7. Adjust PBC for hole
                # If we need to add a hole, we have to scale the system
                # The scaling depends on the type of PBC
                if options["hole"]:
                    if ("square".startswith(options["pbc"]) or
                        "rectangular".startswith(options["pbc"])):
                        scale = 1 + options["hole"] / min(pbcx, pbcy)
                    else:
                        area  = options["hole"]**2/math.cos(math.pi/6)
                        scale = 1+area/(pbcx*pbcy)
                    pbcx, pbcy = scale*pbcx, scale*pbcy

                pbcx = pbcSetX and pbcSetX[0] or pbcx
                pbcy = pbcSetY and pbcSetY[1] or pbcy


                ## 2. Shift of protein relative to the membrane center
                zshift = 0
                if not options["center"]:
                    zshift = -shift[2]
                if options["memshift"]:
                    if options["memshift"] < 0:
                        zshift += options["memshift"] # - max(zip(*prot.coord)[2])
                    else:
                        zshift += options["memshift"] # - min(zip(*prot.coord)[2])

                # Now we center the system in the rectangular
                # brick corresponding to the unit cell
                # If -center is given, also center z in plane
                prot += (0.5*pbcx, 0.5*pbcy, zshift)


            # And we collect the atoms
            protein.atoms.extend(prot.atoms)
            protein.coord.extend(prot.coord)


        # Extract the parts of the protein that are in either leaflet
        prot_up, prot_lo = [], []
        for ix, iy, iz in protein.coord:
            if   iz > 0 and iz <  2.4:
                prot_up.append((ix, iy))
            elif iz < 0 and iz > -2.4:
                prot_lo.append((ix, iy))


        # Current residue ID is set to that of the last atom
        resi = protein.atoms[-1][2]

    atid      = len(protein)+1
    molecules = []

    # The box dimensions are now (likely) set.
    # If a protein was given, it is positioned in the center of the
    # rectangular brick.

    # Set the lattice vectors
    if ("rectangular".startswith(options["pbc"]) or
        "square".startswith(options["pbc"]) or
        "cubic".startswith(options["pbc"])):
        box    = [[pbcx, 0, 0], [0, pbcy, 0], [0, 0, pbcz]]
    elif not lipL:
        # Rhombic dodecahedron with square XY plane
        box    = [[pbcx, 0, 0], [0, pbcy, 0], [0.5*pbcx, 0.5*pbcx, pbcz]]
    elif "hexagonal".startswith(options["pbc"]):
        box    = [[pbcx, 0, 0], [math.sin(math.pi/6)*pbcx, pbcy, 0], [0, 0, pbcz]]
    else: # optimal packing; rhombic dodecahedron with hexagonal XY plane
        box    = [[pbcx, 0, 0], [math.sin(math.pi/6)*pbcx, pbcy, 0], [pbcx/2, pbcy/3, pbcz]]

    # Override lattice vectors if they were set explicitly
    box[0] = pbcSetX or box[0]
    box[1] = pbcSetY or box[1]
    box[2] = pbcSetZ or box[2]

    pbcx, pbcy, pbcz = box[0][0], box[1][1], box[2][2]

    rx, ry, rz = pbcx+1e-8, pbcy+1e-8, pbcz+1e-8


    #################
    ## 2. MEMBRANE ##
    #################

    membrane = Structure()

    if lipL:
        # Lipids are added on grid positions, using the prototypes defined above.
        # If a grid position is already occupied by protein, the position is untagged.

        lipd = lo_lipd

        # Number of lipids in x and y in lower leaflet if there were no solute
        lo_lipids_x = int(pbcx/lipd+0.5)
        lo_lipdx    = pbcx/lo_lipids_x
        lo_rlipx    = range(lo_lipids_x)
        lo_lipids_y = int(pbcy/lipd+0.5)
        lo_lipdy    = pbcy/lo_lipids_y
        lo_rlipy    = range(lo_lipids_y)

        if options["uparea"]:
            lipd = up_lipd

        # Number of lipids in x and y in upper leaflet if there were no solute
        up_lipids_x = int(pbcx/lipd+0.5)
        up_lipdx    = pbcx/up_lipids_x
        up_rlipx    = range(up_lipids_x)
        up_lipids_y = int(pbcy/lipd+0.5)
        up_lipdy    = pbcy/up_lipids_y
        up_rlipy    = range(up_lipids_y)


        # Set up grids to check where to place the lipids
        grid_lo = [[0 for j in lo_rlipy] for i in lo_rlipx]
        grid_up = [[0 for j in up_rlipy] for i in up_rlipx]

        # If there is a protein, mark the corresponding cells as occupied
        if protein:
            # Calculate number density per cell
            for i in prot_lo:
                grid_lo[ int(lo_lipids_x*i[0]/rx)%lo_lipids_x ][ int(lo_lipids_y*i[1]/ry)%lo_lipids_y ] += 1
            for i in prot_up:
                grid_up[ int(up_lipids_x*i[0]/rx)%up_lipids_x ][ int(up_lipids_y*i[1]/ry)%up_lipids_y ] += 1


        # Determine which cells to consider occupied, given the fudge factor
        # The array is changed to boolean type here
        maxd    = float(max([max(i) for i in grid_up+grid_lo]))
        if  maxd == 0:
            if protein:
                print >>sys.stderr, "; The protein seems not to be inside the membrane."
                print >>sys.stderr, "; Run with -orient to put it in."
            maxd = 1


        fudge   = options["fudge"]
        grid_up = [[(j/maxd) <= fudge for j in i] for i in grid_up]
        grid_lo = [[(j/maxd) <= fudge for j in i] for i in grid_lo]


        # If we don't want lipids inside of the protein
        # we also mark everything from the center up to the first cell filled
        if not options["inside"]:

            # Upper leaflet
            marked = [(i, j) for i in up_rlipx for j in up_rlipy if not grid_up[i][j]]
            if marked:
                # Find the center
                cx, cy  = [float(sum(i))/len(marked) for i in zip(*marked)]
                for i, j in marked:
                    md = int(abs(i-cx)+abs(j-cy)) # Manhattan length/distance
                    for f in range(md):
                        ii = int(cx+f*(i-cx)/md)
                        jj = int(cy+f*(j-cy)/md)
                        grid_up[ii][jj] = False

            # Lower leaflet
            marked = [(i, j) for i in lo_rlipx for j in lo_rlipy if not grid_lo[i][j]]
            if marked:
                # Find the center
                cx, cy  = [float(sum(i))/len(marked) for i in zip(*marked)]
                for i, j in marked:
                    md = int(abs(i-cx)+abs(j-cy)) # Manhattan length
                    for f in range(md):
                        ii = int(cx+f*(i-cx)/md)
                        jj = int(cy+f*(j-cy)/md)
                        grid_lo[ii][jj] = False


        # If we make a circular patch, we flag the cells further from the
        # protein or box center than the given radius as occupied.
        if options["disc"]:
            if protein:
                cx, cy = protein.center()[:2]
            else:
                cx, cy = 0.5*pbcx, 0.5*pbcy
            for i in range(len(grid_lo)):
                for j in range(len(grid_lo[i])):
                    if (i*pbcx/lo_lipids_x - cx)**2 + (j*pbcy/lo_lipids_y - cy)**2 > options["disc"]**2:
                        grid_lo[i][j] = False
            for i in range(len(grid_up)):
                for j in range(len(grid_up[i])):
                    if (i*pbcx/up_lipids_x - cx)**2 + (j*pbcy/up_lipids_y - cy)**2 > options["disc"]**2:
                        grid_up[i][j] = False


        # If we need to add a hole, we simply flag the corresponding cells
        # as occupied. The position of the hole depends on the type of PBC,
        # to ensure an optimal arrangement of holes around the protein. If
        # there is no protein, the hole is just put in the center.
        if options["hole"]:
            # Lower leaflet
            if protein:
                if ("square".startswith(options["pbc"]) or
                    "rectangular".startswith(options["pbc"])):
                    hx, hy = (0, 0)
                else:
                    hx, hy = (0, int(lo_lipids_y*math.cos(math.pi/6)/9+0.5))
            else:
                hx, hy = (int(0.5*lo_lipids_x), int(0.5*lo_lipids_y))
            hr = int(options["hole"]/min(lo_lipdx,  lo_lipdy)+0.5)
            ys = int(lo_lipids_x*box[1][0]/box[0][0]+0.5)
            print >>sys.stderr, "; Making a hole with radius %f nm centered at grid cell (%d,%d)"%(options["hole"], hx, hy), hr
            hr -= 1
            for ii in range(hx-hr-1, hx+hr+1):
                for jj in range(hx-hr-1, hx+hr+1):
                    xi, yj = ii, jj
                    if (ii-hx)**2+(jj-hy)**2 < hr**2:
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
            if protein:
                if ("square".startswith(options["pbc"]) or
                    "rectangular".startswith(options["pbc"])):
                    hx, hy = (0, 0)
                else:
                    hx, hy = (0, int(up_lipids_y*math.cos(math.pi/6)/9+0.5))
            else:
                hx, hy = (int(0.5*up_lipids_x), int(0.5*up_lipids_y))
            hr = int(options["hole"]/min(up_lipdx, up_lipdy)+0.5)
            ys = int(up_lipids_x*box[1][0]/box[0][0]+0.5)
            print >>sys.stderr, "; Making a hole with radius %f nm centered at grid cell (%d,%d)"%(options["hole"], hx, hy), hr
            hr -= 1
            for ii in range(hx-hr-1, hx+hr+1):
                for jj in range(hx-hr-1, hx+hr+1):
                    xi, yj = ii, jj
                    if (ii-hx)**2+(jj-hy)**2 < hr**2:
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
        for i in xrange(up_lipids_x):
            for j in xrange(up_lipids_y):
                if grid_up[i][j]:
                    upper.append((random.random(), i*pbcx/up_lipids_x, j*pbcy/up_lipids_y))
        for i in xrange(lo_lipids_x):
            for j in xrange(lo_lipids_y):
                if grid_lo[i][j]:
                    lower.append((random.random(), i*pbcx/lo_lipids_x, j*pbcy/lo_lipids_y))


        # Sort on the random number
        upper.sort()
        lower.sort()


        # Extract coordinates, taking asymmetry in account
        asym  = options["asymmetry"] or 0
        upper = [i[1:] for i in upper[max(0, asym):]]
        lower = [i[1:] for i in lower[max(0, -asym):]]

        print >>sys.stderr, "; X: %.3f (%d bins) Y: %.3f (%d bins) in upper leaflet"%(pbcx, up_lipids_x, pbcy, up_lipids_y)
        print >>sys.stderr, "; X: %.3f (%d bins) Y: %.3f (%d bins) in lower leaflet"%(pbcx, lo_lipids_x, pbcy, lo_lipids_y)
        print >>sys.stderr, "; %d lipids in upper leaflet, %d lipids in lower leaflet"%(len(upper), len(lower))

        # Types of lipids, relative numbers, fractions and numbers

        lipU = lipU or lipL

        # Upper leaflet (+1)
        lipU, numU = zip(*[ parse_mol(i) for i in lipU ])
        totU       = float(sum(numU))
        num_up     = [int(len(upper)*i/totU) for i in numU]
        lip_up     = [l for i, l in zip(num_up, lipU) for j in range(i)]
        leaf_up    = ( 1, zip(lip_up, upper), up_lipd, up_lipdx, up_lipdy)

        # Lower leaflet (-1)
        lipL, numL = zip(*[ parse_mol(i) for i in lipL ])
        totL       = float(sum(numL))
        num_lo     = [int(len(lower)*i/totL) for i in numL]
        lip_lo     = [l for i, l in zip(num_lo, lipL) for j in range(i)]
        leaf_lo    = (-1, zip(lip_lo, lower), lo_lipd, lo_lipdx, lo_lipdy)

        molecules  = zip(lipU, num_up) + zip(lipL, num_lo)

        kick       = options["randkick"]

        # Build the membrane
        for leaflet, leaf_lip, lipd, lipdx, lipdy in [leaf_up, leaf_lo]:
            for lipid, pos in leaf_lip:
                # Increase the residue number by one
                resi += 1
                # Set the random rotation for this lipid
                rangle   = 2*random.random()*math.pi
                rcos     = math.cos(rangle)
                rsin     = math.sin(rangle)
                rcosx    = rcos*lipdx*2/3
                rcosy    = rcos*lipdy*2/3
                rsinx    = rsin*lipdx*2/3
                rsiny    = rsin*lipdy*2/3
                # Fetch the atom list with x, y, z coordinates
                #atoms    = zip(lipidsa[lipid][1].split(), lipidsx[lipidsa[lipid][0]], lipidsy[lipidsa[lipid][0]], lipidsz[lipidsa[lipid][0]])
                # Only keep atoms appropriate for the lipid
                #at, ax, ay, az = zip(*[i for i in atoms if i[0] != "-"])
                at, ax, ay, az = zip(*liplist[lipid].build(diam=lipd))
                # The z-coordinates are spaced at 0.3 nm,
                # starting with the first bead at 0.15 nm
                az       = [ leaflet*(0.5+(i-min(az)))*options["beaddist"] for i in az ]
                xx       = zip( ax, ay )
                nx       = [rcosx*i-rsiny*j+pos[0]+lipdx/2+random.random()*kick for i, j in xx]
                ny       = [rsinx*i+rcosy*j+pos[1]+lipdy/2+random.random()*kick for i, j in xx]
                # Add the atoms to the list
                for i in range(len(at)):
                    atom  = "%5d%-5s%5s%5d"%(resi, lipid, at[i], atid)
                    membrane.coord.append((nx[i], ny[i], az[i]))
                    membrane.atoms.append((at[i], lipid, resi, 0, 0, 0))
                    atid += 1

        # Now move everything to the center of the box before adding solvent
        mz  = pbcz/2
        z   = [ i[2] for i in protein.coord+membrane.coord ]
        mz -= (max(z)+min(z))/2
        protein += (0, 0, mz)
        membrane += (0, 0, mz)


    ################
    ## 3. SOLVENT ##
    ################

    # Charge of the system so far

    last = None
    mcharge = 0
    for j in membrane.atoms:
        if not j[0].strip().startswith('v') and j[1:3] != last:
            mcharge += charges.get(j[1].strip(), 0)
        last = j[1:3]

    last = None
    pcharge = 0
    for j in protein.atoms:
        if not j[0].strip().startswith('v') and j[1:3] != last:
            pcharge += charges.get(j[1].strip(), 0)
        last = j[1:3]

    #mcharge = sum([charges.get(i[0].strip(), 0) for i in set([j[1:3] for j in membrane.atoms])])
    #pcharge = sum([charges.get(i[0].strip(), 0) for i in set([j[1:3] for j in protein.atoms if not j[0].strip().startswith('v')])])



    def _point(y, phi):
        r = math.sqrt(1-y*y)
        return math.cos(phi)*r, y, math.sin(phi)*r


    def pointsOnSphere(n):
        return [_point((2.*k+1)/n-1, k*2.3999632297286531) for k in range(n)]


    if solv:

        # Set up a grid
        d        = 1/options["soldiam"]

        nx, ny, nz = int(1+d*pbcx), int(1+d*pbcy), int(1+d*pbcz)
        dx, dy, dz = pbcx/nx, pbcy/ny, pbcz/nz
        excl, hz  = int(nz*options["solexcl"]/pbcz), int(0.5*nz)

        zshift   = 0
        if membrane:
            memz   = [i[2] for i in membrane.coord]
            midz   = (max(memz)+min(memz))/2
            hz     = int(nz*midz/pbcz)  # Grid layer in which the membrane is located
            zshift = (hz+0.5)*nz - midz # Shift of membrane middle to center of grid layer

        # Initialize a grid of solvent, spanning the whole cell
        # Exclude all cells within specified distance from membrane center
        grid   = [[[i < hz-excl or i > hz+excl for i in xrange(nz)] for j in xrange(ny)] for i in xrange(nx)]

        # Flag all cells occupied by protein or membrane
        for p, q, r in protein.coord+membrane.coord:
            for s, t, u in pointsOnSphere(20):
                x, y, z = p+0.33*s, q+0.33*t, r+0.33*u
                if z >= pbcz:
                    x -= box[2][0]
                    y -= box[2][1]
                    z -= box[2][2]
                if z < 0:
                    x += box[2][0]
                    y += box[2][1]
                    z += box[2][2]
                if y >= pbcy:
                    x -= box[1][0]
                    y -= box[1][1]
                if y < 0:
                    x += box[1][0]
                    y += box[1][1]
                if x >= pbcx:
                    x -= box[0][0]
                if x < 0:
                    x += box[0][0]
                grid[int(nx*x/rx)][int(ny*y/ry)][int(nz*z/rz)] = False

        ##-T grid should be a wrapper around a numpy.ndarray
        # Set the center for each solvent molecule
        kick = options["solrandom"]
        grid = [ (random.random(), (i+0.5+random.random()*kick)*dx, (j+0.5+random.random()*kick)*dy, (k+0.5+random.random()*kick)*dz)
                 for i in xrange(nx) for j in xrange(ny) for k in xrange(nz) if grid[i][j][k] ]

        # Sort on the random number
        grid.sort()

        # 'grid' contains all positions on which a solvent molecule can be placed.
        # The number of positions is taken as the basis for determining the salt concentration.
        # This is fine for simple salt solutions, but may not be optimal for complex mixtures
        # (like when mixing a 1M solution of this with a 1M solution of that

        # First get names and relative numbers for each solvent
        solnames, solnums = zip(*[ parse_mol(i) for i in solv ])
        solnames, solnums = list(solnames), list(solnums)
        totS       = float(sum(solnums))

        # Set the number of ions to add
        nna, ncl = 0, 0
        if options["salt"]:

            # If the concentration is set negative, set the charge to zero
            if options["salt"].startswith("-"):
                charge = 0
                options["salt"] = -float(options["salt"])
            else:
                options["salt"] = float(options["salt"])

            # Determine charge to use, either determined or given on command line
            if options["charge"] != "0":
                charge = (options["charge"] != "auto") and int(options["charge"]) or charge
            else:
                charge = 0

            # Determine number of sodium and chloride to add
            concentration = options["salt"]
            nsol = ("SPC" in solnames and 1 or 4)*len(grid)
            ncl  = max(max(0, charge), int(.5+.5*(concentration*nsol/(27.7+concentration)+charge)))
            nna  = ncl - charge

        # Correct number of grid cells for placement of solvent
        ngrid   = len(grid) - nna - ncl
        num_sol = [int(ngrid*i/totS) for i in solnums]


        # Add salt to solnames and num_sol
        if nna:
            solnames.append("NA+")
            num_sol.append(nna)
            solv.append("NA+")
        if ncl:
            solnames.append("CL-")
            num_sol.append(ncl)
            solv.append("CL-")


        # Names and grid positions for solvent molecules
        solvent    = zip([s for i, s in zip(num_sol, solnames) for j in range(i)], grid)


        # Extend the list of molecules (for the topology)
        molecules.extend(zip(solnames, num_sol))


        # Build the solvent
        sol = []
        for resn, (rndm, x, y, z) in solvent:
            resi += 1
            solmol = solventParticles.get(resn)
            if solmol and len(solmol) > 1:
                # Random rotation (quaternion)
                u,  v,  w       = random.random(), 2*math.pi*random.random(), 2*math.pi*random.random()
                s,  t           = math.sqrt(1-u), math.sqrt(u)
                qw, qx, qy, qz  = s*math.sin(v), s*math.cos(v), t*math.sin(w), t*math.cos(w)
                qq              = qw*qw-qx*qx-qy*qy-qz*qz
                for atnm, (px, py, pz) in solmol:
                    qp = 2*(qx*px + qy*py + qz*pz)
                    rx = x + qp*qx + qq*px + qw*(qy*pz-qz*py)
                    ry = y + qp*qy + qq*py + qw*(qz*px-qx*pz)
                    rz = z + qp*qz + qq*pz + qw*(qx*py-qy*px)
                    sol.append(("%5d%-5s%5s%5d"%(resi%1e5, resn, atnm, atid%1e5), (rx, ry, rz)))
                    atid += 1
            else:
                sol.append(("%5d%-5s%5s%5d"%(resi%1e5, resn, solmol and solmol[0][0] or resn, atid%1e5), (x, y, z)))
                atid += 1
    else:
        solvent, sol = None, []

    return molecules, protein, membrane, solvent, sol, mcharge, pcharge, lipU, lipL, numU, numL, box



def write_all(output, topology, molecules, protein, membrane, solvent, sol, mcharge, pcharge, lipU, lipL, numU, numL, box):
    ## Write the output ##

    grobox = (box[0][0], box[1][1], box[2][2],
              box[0][1], box[0][2], box[1][0],
              box[1][2], box[2][0], box[2][1])

    charge  = mcharge + pcharge
    plen, mlen, slen = 0, 0, 0
    plen = protein and len(protein) or 0
    print >>sys.stderr, "; NDX Solute %d %d" % (1, protein and plen or 0)
    print >>sys.stderr, "; Charge of protein: %f" % pcharge

    mlen = membrane and len(membrane) or 0
    print >>sys.stderr, "; NDX Membrane %d %d" % (1+plen, membrane and plen+mlen or 0)
    print >>sys.stderr, "; Charge of membrane: %f" % mcharge
    print >>sys.stderr, "; Total charge: %f" % charge

    slen = solvent and len(sol) or 0
    print >>sys.stderr, "; NDX Solvent %d %d" % (1+plen+mlen, solvent and plen+mlen+slen or 0)
    print >>sys.stderr, "; NDX System %d %d" % (1, plen+mlen+slen)
    print >>sys.stderr, "; \"I mean, the good stuff is just INSANE\" --Julia Ormond"

    # Open the output stream
    oStream = output and open(output, "w") or sys.stdout

    if output.endswith(".gro"):
        # Print the title
        if membrane.atoms:
            title  = "INSANE! Membrane UpperLeaflet>"+":".join(lipU)+"="+":".join([str(i) for i in numU])
            title += " LowerLeaflet>"+":".join(lipL)+"="+":".join([str(i) for i in numL])

            if protein:
                title = "Protein in " + title
        else:
            title = "Insanely solvated protein."

        print >>oStream, title

        # Print the number of atoms
        print >>oStream, "%5d"%(len(protein)+len(membrane)+len(sol))

        # Print the atoms
        id = 1
        if protein:
            for i in range(len(protein)):
                at, rn, ri = protein.atoms[i][:3]
                x, y, z    = protein.coord[i]
                if rn.endswith('.o'):
                    rn = rn[:-2]
                oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri%1e5, rn, at, id%1e5, x, y, z))
                id += 1
        if membrane:
            for i in range(len(membrane)):
                at, rn, ri = membrane.atoms[i][:3]
                x, y, z    = membrane.coord[i]
                if rn.endswith('.o'):
                    rn = rn[:-2]
                oStream.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(ri%1e5, rn, at, id%1e5, x, y, z))
                id += 1
        if sol:
            # Print the solvent
            print >>oStream, "\n".join([i[0]+"%8.3f%8.3f%8.3f"%i[1] for i in sol])

        # Print the box
        print >>oStream, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f"%grobox
    else:
        # Print the title
        if membrane.atoms:
            title  = "TITLE INSANE! Membrane UpperLeaflet>"+":".join(lipU)+"="+":".join([str(i) for i in numU])
            title += " LowerLeaflet>"+":".join(lipL)+"="+":".join([str(i) for i in numL])
        else:
            title = "TITLE Insanely solvated protein."
        print >>oStream, title

        # Print the box
        print >>oStream, pdbBoxString(box)

        # Print the atoms
        id = 1
        if protein:
            for i in range(len(protein)):
                at, rn, ri = protein.atoms[i][:3]
                x, y, z    = protein.coord[i]
                if rn.endswith('.o'):
                    rn = rn[:-2]
                oStream.write(pdbline%(id%1e5, at, rn, "", ri%1e5, '', 10*x, 10*y, 10*z, 0, 0, ''))
                id += 1
        if membrane:
            for i in range(len(membrane)):
                at, rn, ri = membrane.atoms[i][:3]
                x, y, z    = membrane.coord[i]
                if rn.endswith('.o'):
                    rn = rn[:-2]
                oStream.write(pdbline%(id%1e5, at, rn, "", ri%1e5, '', 10*x, 10*y, 10*z, 0, 0, ''))
                id += 1
        if sol:
            # Print the solvent
            for i in range(len(sol)):
                ri, rn, at, ai = sol[i][0][:5], sol[i][0][5:10], sol[i][0][10:15], sol[i][0][15:20]
                x, y, z    = sol[i][1]
                if rn.endswith('.o'):
                    rn = rn[:-2]
                oStream.write(pdbline%(id%1e5, at.strip(), rn.strip(), "", int(ri)%1e5, '', 10*x, 10*y, 10*z, 0, 0, ''))
                id += 1

    topmolecules = []
    for i in molecules:
        if i[0].endswith('.o'):
            topmolecules.append(tuple([i[0][:-2]]+list(i[1:])))
        else:
            topmolecules.append(i)

    if topology:
        # Write a rudimentary topology file
        with open(topology, "w") as top:
            print >>top, '#include "martini.itp"\n'
            print >>top, '[ system ]\n; name\n%s\n\n[ molecules ]\n; name  number'%title
            if protein:
                print >>top, "%-10s %5d"%("Protein", 1)
            print >>top, "\n".join("%-10s %7d"%i for i in topmolecules)
    else:
        print >>sys.stderr, "\n".join("%-10s %7d"%i for i in topmolecules)


def insane(**options):

    ## PROTEINS
    # Collect and build input structures (proteins)
    # Determine protein/hole layout
    # Collect and build lipid/solvent dictionaries

    ## PBC
    # Determine size of system

    # MEMBRANE
    # Make grid - two leaflets
    # Mask grid occupied cells
    # Determine lipid compositions per leaflet
    # Build membrane

    ## SOLVENT
    # Determine solvent composition
    # Build solvent

    return 0
