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

from __future__ import print_function

import os
import sys
import random
import collections
import numpy as np

from . import linalg
from . import lipids
from .converters import *
from .constants import d2r
from ._data import SOLVENTS, CHARGES, APOLARS

try:
    range = xrange
except NameError:
    pass

RTOL = 1e-8

# Set the random seed.
# The seed is set to an arbitary value set in the INSANE_SEED environment
# variable. If the environment variable is not set, then the system time is
# used to set the seed.
random.seed(os.environ.get('INSANE_SEED', None))


class PBCException(Exception):
    pass

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
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1"
pdbline = "ATOM  %5i  %-3s %4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s  "


def pdbBoxString(box):
    # Box vectors
    u, v, w  = box

    # Box vector lengths
    nu, nv, nw = [np.sqrt(linalg.norm2(i)) for i in (u, v, w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or np.arccos(linalg.cos_angle(v, w))/d2r
    beta  = nu*nw == 0 and 90 or np.arccos(linalg.cos_angle(u, w))/d2r
    gamma = nu*nv == 0 and 90 or np.arccos(linalg.cos_angle(u, v))/d2r

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


class PBC(object):
    def __init__(self, shape=None, box=None, xyz=None, distance=None, membrane=None, protein=None, disc=None, hole=None):
        self.box = None

        if box:
            #print("Setting box from given box definition. Neglecting all other PBC options!")
            self.box = np.array(box).reshape((3,3))
            return

        x, y, z = None, None, None
        if xyz and any(xyz):
            x, y, z = xyz
            if type(x) in (float, int):
                x = (x, 0, 0)
            if type(y) in (float, int):
                y = (0, y, 0)
            if type(z) in (float, int):
                z = (0, 0, z)
            if all(xyz):
                #print("Setting box from x/y/z vectors provided. Neglecting all other PBC options!")
                self.box = np.array([x,y,z])
                return
        xyz = (x, y, z)

        # Not having a distance here is an error
        if distance is None:
            raise PBCException("Cannot set PBC if size of system is not specified.")

        # Having a distance at zero without a solute for the extent is also an error
        if not distance and not (protein or hole or disc):
            raise PBCException("Having a PBC distance of 0 without having a "
                               "solute/hole/disc results in a box with zero volume.")

        if type(distance) in (float, int):
            scale = zscale = distance
        else:
            scale, zscale = distance

        if disc is not None:
            scale += 2*disc
        elif hole is not None:
            scale += 2*hole

        if protein:
            # Assume proteins are centered - use min/max
            pmin, pmax = protein[0].fun(min), protein[0].fun(max)
            prng       = (pmax[0]-pmin[0], pmax[1]-pmin[1], pmax[2]-pmin[2])
            zmax = max([ p.fun(max)[2] for p in protein ])
            zmin = min([ p.fun(min)[2] for p in protein ])
            zscale += zmax - zmin

        if x and y and not z:
            #print("Setting box from x/y vectors provided and distance set over z "
            #      "(including solute dimensions). Neglecting other PBC options!")
            self.box = np.array([x,y,[0,0,zscale]])
            return

        # After this point, we need to finalize the box by (re)setting
        # x, y and/or z, according to the xyz option.

        if "cubic".startswith(shape) or ("square".startswith(shape) and not membrane):
            if protein:
                scale += protein[0].diam()
            self.box = np.eye(3) * scale
        elif "rectangular".startswith(shape):
            if protein:
                self.box = np.diag(prng) + scale
            elif disc or hole:
                self.box = np.array([[scale, 0, 0],
                                     [0, scale, 0],
                                     [0, 0, zscale]])
            else:
                self.box = np.diag((scale,scale,zscale))
        elif not membrane:
            if protein:
                scale += protein[0].diam()
            # Replace with hexagonal prism with long axis along X
            # -- this one makes no sense at all!
            if "hexagonal".startswith(shape) and not protein:
                self.box = np.array([[1, 0, 0],
                                     [0, np.sqrt(0.75), 0],
                                     [0.5, 0.5, 1]])*scale
                self.box[2,2] = zscale
            else:
            # No membrane, no cubic/rectangular box, returning a rhombic dodecahedron
            # (matches for "dodecahedron", "optimal")
#            self.box = np.array([[1, 0, 0],
#                                 [0.5, np.sqrt(0.75), 0],
#                                 [0.5, np.sqrt(3)/6, np.sqrt(6)/3]])*scale
                self.box = np.array([[1, 0, 0],
                                     [0, 1, 0],
                                     [0.5, 0.5, np.sqrt(0.5)]])*scale
                #self.box[2,2] = zscale
        # Only the membrane cases left here
        else:
            if protein and not (disc or hole):
                #scale  += protein[0].diamxy()
                scale  += prng[0]
            if "square".startswith(shape):
                self.box = np.array([[scale, 0, 0],
                                     [0, scale, 0],
                                     [0, 0, zscale]])
            elif "optimal".startswith(shape):
                #print("Warning: this box may be too skewed for Gromacs.")
                #self.box = np.array([[scale, 0, 0],
                #                     [scale/2, scale*np.sqrt(0.75), 0],
                #                     [scale/2, scale*np.sqrt(3)/6, 0]])
                #self.box[2,2] = np.sqrt(zscale**2 - (self.box[2,:]**2).sum())
                self.box = np.array([[scale, 0, 0],
                                     [0, scale, 0],
                                     [scale/2, scale/2, 0]])
                self.box[2,2] = np.sqrt(zscale**2 - (self.box[2,:]**2).sum())
            else:
                self.box = np.array([[scale, 0, 0],
                                     [scale/2, scale*np.sqrt(0.75), 0],
                                     [0, 0, zscale]])

        for i, v in enumerate(xyz):
            if v:
                self.box[i,:] = v

        self.box = self.box.astype(np.float64)
        return

    @property
    def x(self):
        return self.box[0, 0]

    @x.setter
    def x(self, value):
        self.box[0, 0] = value

    @property
    def y(self):
        return self.box[1, 1]

    @y.setter
    def y(self, value):
        self.box[1, 1] = value

    @property
    def z(self):
        return self.box[2,2]

    @z.setter
    def z(self, value):
        self.box[2, 2] = value

    @property
    def rx(self):
        return self.box[0,0] + RTOL

    @property
    def ry(self):
        return self.box[1,1] + RTOL

    @property
    def rz(self):
        return self.box[2,2] + RTOL


class Structure(object):
    def __init__(self, filename=None):
        self.title   = ""
        self.atoms   = []
        self._coord  = None
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


    def __nonzero__(self):
        return bool(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __iadd__(self, s):
        if self.coord.shape[0]:
            self.coord += s ###
        return self

    def __add__(self, other):
        if hasattr(other, 'atoms') and hasattr(other, 'coord'):
            result = self.__class__()
            result.atoms.extend(self.atoms)
            result.atoms.extend(other.atoms)
            result.coord = np.concatenate((
                self.coord.reshape((-1,3)),
                other.coord.reshape((-1,3)))) ###
            return result
        raise TypeError('Cannot add {} to {}'
                        .format(self.__class__, other.__class__))

    def __iter__(self):
        atom_enumeration = enumerate(zip(self.atoms, self.coord), start=1)
        for idx, (atom, (x, y, z)) in atom_enumeration:
            atname, resname, resid = atom[:3]
            if resname.endswith('.o'):
                resname = rn[:-2]
            yield idx, atname, resname, resid, x, y, z

    @property
    def coord(self):
        if self._coord is None:
            self._coord = np.array([i[4:7] for i in self.atoms]).reshape((-1,3))
        return self._coord

    @coord.setter
    def coord(self, other):
        self._coord = np.array(other).reshape((-1,3))

    @property
    def charge(self):
        last = None
        charge = 0
        for j in self.atoms:
            if not j[0].strip().startswith('v') and j[1:3] != last:
                charge += CHARGES.get(j[1].strip(), 0)
            last = j[1:3]
        return charge

    @property
    def center(self):
        if self._center is None:
            self._center = self.coord.mean(axis=0)
        return self._center

    @center.setter
    def center(self, other):
        s = other - self.coord.mean(axis=0)
        self.coord += s ###

    def diam(self):
        if np.any(self._center):
            self.center = (0, 0, 0)
        return 2*np.sqrt(max([i*i+j*j+k*k for i, j, k in self.coord]))

    def diamxy(self):
        if np.any(self._center):
            self.center = (0, 0, 0)
        return 2*np.sqrt(max([i*i+j*j for i, j, k in self.coord]))

    def areaxy(self, lowerbound=-np.inf, upperbound=np.inf, spacing=0.1):
        mask = (self.coord[:,2] > lowerbound) & (self.coord[:,2] < upperbound)
        points = self.coord[mask, :2]
        # The magic number factor 1.1 is not critical at all
        # Just a number to set a margin to the bounding box and 
        # have all points fall within the boundaries
        bbmin, bbmax = 1.1*points.min(axis=0), 1.1*points.max(axis=0)
        size = bbmax - bbmin
        cells = (size / spacing + 0.5).astype('int')
        # Grid points over bounding box with specified spacing
        grid = np.mgrid[bbmin[0]:bbmax[0]:(cells[0]*1j),
                        bbmin[1]:bbmax[1]:(cells[1]*1j)].reshape((2,-1)).T
        # Occupied cells is approximately equal to grid points within
        # gridspacing distance of points
        occupied = occupancy(grid, points, spacing)
        # The occupied area follows from the fraction of occupied
        # cells times the area spanned by the bounding box
        return size[0]*size[1]*sum(occupied > 0)/occupied.size


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
                phobic[jx][jy][jz] += (i[1].strip() in APOLARS)

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
        self.center = (-sxm, -sym, -szm)
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
        self.coord = np.array([(ux*i+uy*j+uz*k, vx*i+vy*j+vz*k, wx*i+wy*j+wz*k)
                               for i, j, k in self.coord])


    def rotate(self, what):
            if what == "princ":
                self.rotate_princ()
            ## ii. Randomly
            elif what == "random":
                self.rotate_random()
            ## iii. Specifically
            elif what:
                self.rotate_degrees(float(what))


    def rotate_princ(self):
        R = np.linalg.eig(np.dot(self.coord[:,:2].T,self.coord[:,:2]))
        self.coord[:,:2] = np.dot(self.coord[:,:2], R[1][:,np.argsort(R[0])[::-1]])
        return


    def rotate_random(self):
        ux   = np.cos(random.random()*2*np.pi)
        uy   = np.sqrt(1-ux*ux)
        self.coord[:,:2] = np.dot(self.coord[:,:2],[[ux,-uy],[uy,ux]])


    def rotate_degrees(self, angle):
        ux   = np.cos(angle*np.pi/180.)
        uy   = np.sin(angle*np.pi/180.)
        self.coord[:,:2] = np.dot(self.coord[:,:2], [[ux, -uy],[uy, ux]])


def _point(y, phi):
    r = np.sqrt(1-y*y)
    return np.cos(phi)*r, y, np.sin(phi)*r


def pointsOnSphere(n):
    return np.array([_point((2.*k+1)/n-1, k*2.3999632297286531) for k in range(n)])


# Mean of deviations from initial value
def meand(v):
    return sum([i-v[0] for i in v])/len(v)


# Sum of squares/crossproducts of deviations
def ssd(u, v):
    return sum([(i-u[0])*(j-v[0]) for i, j in zip(u, v)])/(len(u)-1)


def occupancy(grid, points, spacing=0.01):
    """Return a vector with the occupancy of each grid point for 
    given array of points"""
    distances = ((grid[:,None,:] - points[None,:,:])**2).sum(axis=2)
    occupied = (distances < spacing).sum(axis=1)
    return occupied


def determine_molecule_numbers(total, molecules, absolute, relative):
    """Determine molecule numbers for given total, 
    absolute and relative numbers""" 
    weight = sum(relative)
    if not any(absolute):
        # Only relative numbers
        numbers = [int(total*i/weight) for i in relative]
    elif any(relative):
        # Absolute numbers and fill the rest with relative numbers
        rest = total - sum(absolute)
        numbers = [int(rest*i/weight) if i else j 
                   for i,j in zip(relative, absolute)]
    else:
        # Only absolute numbers
        numbers = absolute
    return list(zip(molecules, numbers))


def resize_pbc_for_lipids(pbc, relL, relU, absL, absU,
                          uparea, area, hole, proteins):
    """
    Adapt the size of the box to accomodate the lipids.

    The PBC is changed **on place**.
    """
    if any(relL) and any(relU):
        # Determine box from size
        # If one leaflet is defined with an absolute number of lipids
        # then the other leaflet (with relative numbers) will follow
        # from that.
        # box/d/x/y/z needed to define unit cell

        # box is set up already...
        # yet it may be set to 0, then there is nothing we can do.
        if 0 in (pbc.x, pbc.y, pbc.z):
            raise PBCException('Not enough information to set the box size.')
    elif any(absL) or any(absU):
        # All numbers are absolute.. determine size from number of lipids
        # Box x/y will be set, d/dz/z is needed to set third box vector.
        # The area is needed to determine the size of the x/y plane OR
        # the area will be SET if a box is given.

        if pbc.z == 0:
            raise PBCException('Not enough information to set the box size.')

        if 0 in (pbc.x, pbc.y):
            # We do not know what size the box should be. Let's X and Y should
            # be the same.
            pbc.x = pbc.y = 1
        # A scaling factor is needed for the box
        # This is the area for the given number of lipids
        upsize = sum(absU) * uparea
        losize = sum(absL) * area
        # This is the size of the hole, going through both leaflets
        holesize = np.pi * hole ** 2
        # This is the area of the PBC xy plane
        xysize = pbc.x * pbc.y
        # This is the total area of the proteins per leaflet (IMPLEMENT!)
        psize_up = sum([p.areaxy(0, 2.4) for p in proteins])
        psize_lo = sum([p.areaxy(-2.4, 0) for p in proteins])
        # So this is unavailable:
        unavail_up = holesize + psize_up
        unavail_lo = holesize + psize_lo
        # This is the current area marked for lipids
        # xysize_up = xysize - unavail_up
        # xysize_lo = xysize - unavail_lo
        # This is how much the unit cell xy needs to be scaled
        # to accomodate the fixed amount of lipids with given area.
        upscale = (upsize + unavail_up)/xysize
        loscale = (losize + unavail_lo)/xysize
        area_scale = max(upscale, loscale)
        aspect_ratio = pbc.x / pbc.y
        scale_x = np.sqrt(area_scale / aspect_ratio)
        scale_y = np.sqrt(area_scale / aspect_ratio)
        pbc.box[:2,:] *= math.sqrt(area_scale)


def setup_solvent(pbc, protein, membrane, options):
    # Charge of the system so far

    if not options["solvent"]:
        return Structure(), []

    solv = options["solvent"]
    charge = membrane.charge + protein.charge

    # Set up a grid
    d        = 1/options["soldiam"]

    nx, ny, nz = int(1+d*pbc.x), int(1+d*pbc.y), int(1+d*pbc.z)
    dx, dy, dz = pbc.x/nx, pbc.y/ny, pbc.z/nz
    excl, hz  = int(nz*options["solexcl"]/pbc.z), int(0.5*nz)

    zshift   = 0
    if membrane:
        memz   = [i[2] for i in membrane.coord]
        midz   = (max(memz)+min(memz))/2
        hz     = int(nz*midz/pbc.z)  # Grid layer in which the membrane is located
        zshift = (hz+0.5)*nz - midz # Shift of membrane middle to center of grid layer

    # Initialize a grid of solvent, spanning the whole cell
    # Exclude all cells within specified distance from membrane center
    grid   = [[[i < hz-excl or i > hz+excl for i in range(nz)] for j in range(ny)] for i in range(nx)]

    # Flag all cells occupied by protein or membrane
    for coord in (protein+membrane).coord:
        for x, y, z in coord + 0.33*pointsOnSphere(20):
            if z >= pbc.z:
                x -= pbc.box[2,0]
                y -= pbc.box[2,1]
                z -= pbc.box[2,2]
            if z < 0:
                x += pbc.box[2,0]
                y += pbc.box[2,1]
                z += pbc.box[2,2]
            if y >= pbc.y:
                x -= pbc.box[1,0]
                y -= pbc.box[1,1]
            if y < 0:
                x += pbc.box[1,0]
                y += pbc.box[1,1]
            if x >= pbc.x:
                x -= pbc.box[0,0]
            if x < 0:
                x += pbc.box[0,0]
            grid[int(nx*x/pbc.rx)][int(ny*y/pbc.ry)][int(nz*z/pbc.rz)] = False

    ##-T grid should be a wrapper around a numpy.ndarray
    # Set the center for each solvent molecule
    kick = options["solrandom"]
    grid = [ (random.random(), (i+0.5+random.random()*kick)*dx, (j+0.5+random.random()*kick)*dy, (k+0.5+random.random()*kick)*dz)
             for i in range(nx) for j in range(ny) for k in range(nz) if grid[i][j][k] ]

    # Sort on the random number
    grid.sort()

    # 'grid' contains all positions on which a solvent molecule can be placed.
    # The number of positions is taken as the basis for determining the salt concentration.
    # This is fine for simple salt solutions, but may not be optimal for complex mixtures
    # (like when mixing a 1M solution of this with a 1M solution of that

    # First get names and relative numbers for each solvent
    solnames, solabs, solnums = zip(*solv)
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
    solvent    = list(zip([s for i, s in zip(num_sol, solnames) for j in range(i)], grid))

    # Build the solvent
    resi = 0
    if protein:
        resi = protein.atoms[-1][2]
    if membrane:
        resi = membrane.atoms[-1][2]
    sol = Structure()
    solcoord = []
    for resn, (rndm, x, y, z) in solvent:
        resi += 1
        solmol = SOLVENTS.get(resn)
        if solmol and len(solmol) > 1:
            # Random rotation (quaternion)
            u,  v,  w       = random.random(), 2*np.pi*random.random(), 2*np.pi*random.random()
            s,  t           = np.sqrt(1-u), np.sqrt(u)
            qw, qx, qy, qz  = s*np.sin(v), s*np.cos(v), t*np.sin(w), t*np.cos(w)
            qq              = qw*qw-qx*qx-qy*qy-qz*qz
            for atnm, (px, py, pz) in solmol:
                qp = 2*(qx*px + qy*py + qz*pz)
                rx = x + qp*qx + qq*px + qw*(qy*pz-qz*py)
                ry = y + qp*qy + qq*py + qw*(qz*px-qx*pz)
                rz = z + qp*qz + qq*pz + qw*(qx*py-qy*px)
                sol.atoms.append((atnm, resn, resi, 0, 0, 0))
                solcoord.append((rx, ry, rz))
        else:
            sol.atoms.append((solmol and solmol[0][0] or resn,
                              resn, resi,
                              0, 0, 0))
            solcoord.append((x, y, z))
    sol.coord = solcoord

    return sol, zip(solnames, num_sol)


def setup_membrane(pbc, protein, lipid, options):
    membrane = Structure()
    molecules = []

    lower, upper = lipid
    lipL, absL, relL = lower
    lipU, absU, relU = upper

    if not any((absL, relL, absU, relU)):
        return membrane, molecules

    lo_lipd  = np.sqrt(options["area"])
    up_lipd = np.sqrt(options["uparea"])

    num_up, num_lo = [], []

    # Lipids are added on grid positions, using the prototypes defined above.
    # If a grid position is already occupied by protein, the position is untagged.
    
    lipd = lo_lipd

    # Number of lipids in x and y in lower leaflet if there were no solute
    lo_lipids_x = int(pbc.x/lo_lipd+0.5)
    lo_lipdx    = pbc.x/lo_lipids_x
    lo_rlipx    = list(range(lo_lipids_x))
    lo_lipids_y = int(pbc.y/lo_lipd+0.5)
    lo_lipdy    = pbc.y/lo_lipids_y
    lo_rlipy    = list(range(lo_lipids_y))
    
    if options["uparea"]:
        lipd = up_lipd

    # Number of lipids in x and y in upper leaflet if there were no solute
    up_lipids_x = int(pbc.x/up_lipd+0.5)
    up_lipdx    = pbc.x/up_lipids_x
    up_rlipx    = list(range(up_lipids_x))
    up_lipids_y = int(pbc.y/up_lipd+0.5)
    up_lipdy    = pbc.y/up_lipids_y
    up_rlipy    = list(range(up_lipids_y))

    # Set up grids to check where to place the lipids
    grid_lo = [[0 for j in lo_rlipy] for i in lo_rlipx]
    grid_up = [[0 for j in up_rlipy] for i in up_rlipx]

    ## NEW GRID
    lo_gx = slice(0, pbc.x, int(pbc.x/lo_lipd + 0.5)*1j)
    lo_gy = slice(0, pbc.x, int(pbc.x/lo_lipd + 0.5)*1j)
    up_gx = slice(0, pbc.x, int(pbc.x/up_lipd + 0.5)*1j)
    up_gy = slice(0, pbc.x, int(pbc.x/up_lipd + 0.5)*1j)
    grid_l = np.mgrid[lo_gx, lo_gy][:,:-1,:-1].reshape((2,-1)).T
    grid_u = np.mgrid[up_gx, up_gy][:,:-1,:-1].reshape((2,-1)).T
    # If there is a protein, mark the corresponding POINTS as occupied
    occupied_lo = np.zeros(grid_l.shape[0])
    occupied_up = np.zeros(grid_u.shape[0])
    if protein:
        upmask = (protein.coord[:,2] <  2.4) & (protein.coord[:,2] > 0)
        lomask = (protein.coord[:,2] > -2.4) & (protein.coord[:,2] < 0)
        occupied_lo += occupancy(grid_l, protein.coord[lomask,:2], lo_lipd)
        occupied_up += occupancy(grid_u, protein.coord[upmask,:2], up_lipd)
        maxd = max(occupied_lo.max(), occupied_up.max())
        if maxd:
            occupied_up = (occupied_up/maxd) > options["fudge"]
            occupied_lo = (occupied_lo/maxd) > options["fudge"]
    #print(absU, sum(~occupied_up.astype('bool')), 
    #      absL, sum(~occupied_lo.astype('bool')))

    ## OLD GRID
        
    maxd = 1

    # If there is a protein, mark the corresponding cells as occupied
    if protein:
        # Extract the parts of the protein that are in either leaflet
        mem_mask_up = (0 < protein.coord[:,2]) & (protein.coord[:,2] < 2.4)
        mem_mask_lo = (0 > protein.coord[:,2]) & (protein.coord[:,2] > -2.4)
        prot_up = protein.coord[mem_mask_up, :2]
        prot_lo = protein.coord[mem_mask_lo, :2]

        # Calculate number density per cell
        for i in prot_lo:
            grid_lo[ int(lo_lipids_x*i[0]/pbc.rx)%lo_lipids_x ][ int(lo_lipids_y*i[1]/pbc.ry)%lo_lipids_y ] += 1
        for i in prot_up:
            grid_up[ int(up_lipids_x*i[0]/pbc.rx)%up_lipids_x ][ int(up_lipids_y*i[1]/pbc.ry)%up_lipids_y ] += 1

        # Determine which cells to consider occupied, given the fudge factor
        # The array is changed to boolean type here
        maxd = float(max([max(i) for i in grid_up+grid_lo]))
        if  maxd == 0:
            print("; The protein seems not to be inside the membrane.", file=sys.stderr)
            print("; Run with -orient to put it in.", file=sys.stderr)
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
            cx, cy = protein.center[:2]
        else:
            cx, cy = 0.5*pbc.x, 0.5*pbc.y
        for i in range(len(grid_lo)):
            for j in range(len(grid_lo[i])):
                if (i*pbc.x/lo_lipids_x - cx)**2 + (j*pbc.y/lo_lipids_y - cy)**2 > options["disc"]**2:
                    grid_lo[i][j] = False
        for i in range(len(grid_up)):
            for j in range(len(grid_up[i])):
                if (i*pbc.x/up_lipids_x - cx)**2 + (j*pbc.y/up_lipids_y - cy)**2 > options["disc"]**2:
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
                hx, hy = (0, int(lo_lipids_y*np.cos(np.pi/6)/9+0.5))
        else:
            hx, hy = (int(0.5*lo_lipids_x), int(0.5*lo_lipids_y))
        hr = int(options["hole"]/min(lo_lipdx,  lo_lipdy)+0.5)
        ys = int(lo_lipids_x*pbc.box[1,0]/pbc.box[0,0]+0.5)
        print("; Making a hole with radius %f nm centered at grid cell (%d,%d)"%(options["hole"], hx, hy), hr, file=sys.stderr)
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
                hx, hy = (0, int(up_lipids_y*np.cos(np.pi/6)/9+0.5))
        else:
            hx, hy = (int(0.5*up_lipids_x), int(0.5*up_lipids_y))
        hr = int(options["hole"]/min(up_lipdx, up_lipdy)+0.5)
        ys = int(up_lipids_x*pbc.box[1,0]/pbc.box[0,0]+0.5)
        print("; Making a hole with radius %f nm centered at grid cell (%d,%d)"%(options["hole"], hx, hy), hr, file=sys.stderr)
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
    for i in range(up_lipids_x):
        for j in range(up_lipids_y):
            if grid_up[i][j]:
                upper.append((random.random(), i*pbc.x/up_lipids_x, j*pbc.y/up_lipids_y))
    for i in range(lo_lipids_x):
        for j in range(lo_lipids_y):
            if grid_lo[i][j]:
                lower.append((random.random(), i*pbc.x/lo_lipids_x, j*pbc.y/lo_lipids_y))

    # Sort on the random number
    upper.sort()
    lower.sort()

    # Extract coordinates, taking asymmetry in account
    asym  = options["asymmetry"] or 0
    upper = [i[1:] for i in upper[max(0, asym):]]
    lower = [i[1:] for i in lower[max(0, -asym):]]

    print("; X: %.3f (%d bins) Y: %.3f (%d bins) in upper leaflet"%(pbc.x, up_lipids_x, pbc.y, up_lipids_y), file=sys.stderr)
    print("; X: %.3f (%d bins) Y: %.3f (%d bins) in lower leaflet"%(pbc.x, lo_lipids_x, pbc.y, lo_lipids_y), file=sys.stderr)
    print("; %d lipids in upper leaflet, %d lipids in lower leaflet"%(len(upper), len(lower)), file=sys.stderr)

    ##> Types of lipids, relative numbers, fractions and numbers

    # Upper leaflet (+1)
    numbers = determine_molecule_numbers(len(upper), lipU, absU, relU)
    lip_up     = [l for l, i in numbers for j in range(i)]
    leaf_up    = ( 1, zip(lip_up, upper), up_lipd, up_lipdx, up_lipdy)
    molecules.extend(numbers)

    # Lower leaflet (-1)
    numbers = determine_molecule_numbers(len(lower), lipL, absL, relL)
    lip_lo = [l for l, i in numbers for j in range(i)]
    leaf_lo = (-1, zip(lip_lo, lower), lo_lipd, lo_lipdx, lo_lipdy)
    molecules.extend(numbers)

    ##< Done determining numbers per lipid

    ##> Building lipids

    kick       = options["randkick"]

    mematoms = []
    memcoords = []

    # Build the membrane

    ## ==> LIPID  BOOKKEEPING:
    # Read lipids defined in insane
    liplist = lipids.get_lipids()
    # Then add lipids from file
    liplist.add_from_files(options["molfile"])
    # Last, add lipids from command line
    liplist.add_from_def(options["lipnames"], options["lipheads"], options["liplinks"],
                         options["liptails"], options["lipcharge"])

    if protein:
        resi = protein.atoms[-1][2]
    else:
        resi = 0

    for leaflet, leaf_lip, lipd, lipdx, lipdy in [leaf_up, leaf_lo]:
        for lipid, pos in leaf_lip:
            # Increase the residue number by one
            resi += 1
            # Set the random rotation for this lipid
            rangle   = 2*random.random()*np.pi
            rcos     = np.cos(rangle)
            rsin     = np.sin(rangle)
            rcosx    = rcos*lipdx*2/3
            rcosy    = rcos*lipdy*2/3
            rsinx    = rsin*lipdx*2/3
            rsiny    = rsin*lipdy*2/3
            # Fetch the atom list with x, y, z coordinates
            at, ax, ay, az = zip(*liplist[lipid].build(diam=lipd))
            # The z-coordinates are spaced at 0.3 nm,
            # starting with the first bead at 0.15 nm
            az = [ leaflet*(0.5+(i-min(az)))*options["beaddist"] for i in az ]
            xx = np.array((ax, ay)).T
            nx = np.dot(xx,(rcosx, -rsiny)) + pos[0] + lipdx/2 + [ kick*random.random() for i in az ]
            ny = np.dot(xx,(rsinx, rcosy)) + pos[1] + lipdy/2 + [ kick*random.random() for i in az ]
            # Add the atoms to the list
            memcoords.extend([(nx[i], ny[i], az[i]) for i in range(len(at))])
            mematoms.extend([(at[i], lipid, resi, 0, 0, 0) for i in range(len(at))])

    ##< Done building lipids

    membrane.coord = memcoords
    membrane.atoms = mematoms

    return membrane, molecules


def old_main(**options):

    molecules = []

    #########################################
    ## I. PROTEIN and other macromolecules ##
    #########################################

    # Read in the structures (if any)
    tm = [ Structure(i) for i in options["solute"] ]

    xshifts  = [0] # Shift in x direction per protein

    if tm:
        molecules.append(('Protein', len(tm)))

    for prot in tm:
        # Center the protein and store the shift
        shift = prot.center
        prot.center = (0, 0, 0)

        ## 1. Orient with respect to membrane
        # Orient the protein according to the TM region, if requested
        # This doesn't actually work very well...
        if options["orient"]:
            prot.orient(options["origriddist"], options["oripower"])

        ## 4. Orient the protein in the xy-plane
        ## i. According to principal axes and unit cell
        prot.rotate(options["rotate"])

        ## 5. Determine the minimum and maximum x and y of the protein
        pmin, pmax = prot.coord.min(axis=0), prot.coord.max(axis=0)

        # At this point we should shift the subsequent proteins such
        # that they end up at the specified distance, in case we have
        # a number of them to do
        # y-shift is always -ycenter
        # x-shift is -xmin+distance+xmax(current)
        xshifts.append(xshifts[-1]+pmax[0]+(options["distance"] or 0))

        ## 2. Shift of protein relative to the membrane center
        zshift = options["memshift"]
        if not options["center"]:
            zshift -= shift[2]

        # Now we center the system in the rectangular
        # brick corresponding to the unit cell
        # If -center is given, also center z in plane
        prot += (0, 0, zshift)

        # The z position is now correct with respect to the membrane
        # at z = 0. The x/y need to be set still

    #############
    ## II. PBC ##
    #############

    # Periodic boundary conditions
    if options["pbc"] == 'keep' and tm:
        box = tm[0].box
    else:
        box = options.get("box")

    zdist = options["zdistance"]
    if zdist == None:
        zdist = options["distance"]

    # Set up base PBC
    # Override where needed to accomodate additional components
    # box/shape are final - if these are given and a solute does
    # not fit in it raises an exception
    pbc = PBC(shape=options["pbc"], box=box,
              distance=(options["distance"], zdist),
              xyz=(options["xvector"], options["yvector"], options["zvector"]),
              disc=options["disc"], hole=options["hole"],
              membrane=options["lower"], protein=tm)


    #################
    ## III. LIPIDS ##
    #################

    lipL      = options["lower"]
    lipU      = options["upper"]
    relU, relL, absU, absL = [], [], [], []
    if lipL:
        lipU = lipU or lipL
        lipL, absL, relL = zip(*lipL)
        totL       = float(sum(relL))
        lipU, absU, relU = zip(*lipU)
        totU       = float(sum(relU))
    else:
        options["solexcl"] = -1

    if options["uparea"] is None:
        options["uparea"] = options["area"]

    resize_pbc_for_lipids(pbc=pbc, relL=relL, relU=relU, absL=absL, absU=absU,
                          uparea=options["uparea"], area=options["area"],
                          hole=options["hole"], proteins=tm)

    ##################
    ## IV. MEMBRANE ##
    ##################

    ## 1. Proteins

    # Now that PBC is set, we can shift the proteins
    for xshft, prot in zip(xshifts,tm):
        # Half the distance should be added to xshft
        # to center the whole lot.
        #xshft += options["distance"]/2
        xshft = pbc.x/2
        prot.coord += (xshft, pbc.y/2, (not lipL)*pbc.z/2)

    prot_coord = []
    protein  = Structure()
    for prot in tm:
        # And we collect the atoms
        protein.atoms.extend(prot.atoms)
        prot_coord.append(prot.coord)
    if tm:
        protein.coord = np.concatenate(prot_coord)

    # Current residue ID is set to that of the last atom
    resi = 0 
    if protein:
        resi = protein.atoms[-1][2]

    ## 2. Lipids

    lipid = ((lipL, absL, relL), (lipU, absU, relU))
    membrane, added = setup_membrane(pbc, protein, lipid, options)
    molecules.extend(added)

    if added:
        # Now move everything to the center of the box before adding solvent
        mz  = pbc.z/2
        z   = (protein+membrane).coord[:,2]
        mz -= (max(z)+min(z))/2
        protein.coord += (0, 0, mz)
        membrane.coord += (0, 0, mz)

    if membrane:
        resi = membrane.atoms[-1][2]

    ################
    ## 3. SOLVENT ##
    ################


    solvent, added = setup_solvent(pbc, protein, membrane, options)
    molecules.extend(added)

    return (molecules, protein, membrane, solvent, lipid, pbc.box)



def write_gro(outfile, title, atoms, box):
    """
    Write a GRO file.

    Parameters
    ----------
    outfile
        The stream to write in.
    title
        The title of the GRO file. Must be a single line.
    atoms
        An instance of Structure containing the atoms to write.
    box
        The periodic box as a 3x3 matrix.
    """
    # Print the title
    print(title, file=outfile)

    # Print the number of atoms
    print("{:5d}".format(len(atoms)), file=outfile)

    # Print the atoms
    atom_template = "{:5d}{:<5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}"
    for idx, atname, resname, resid, x, y, z in atoms:
        print(atom_template
              .format(int(resid % 1e5), resname, atname, int(idx % 1e5),
                      x, y, z),
              file=outfile)

    # Print the box
    grobox = (box[0][0], box[1][1], box[2][2],
              box[0][1], box[0][2], box[1][0],
              box[1][2], box[2][0], box[2][1])
    box_template = '{:10.5f}' * 9
    print(box_template.format(*grobox), file=outfile)


def write_pdb(outfile, title, atoms, box):
    """
    Write a PDB file.

    Parameters
    ----------
    outfile
        The stream to write in.
    title
        The title of the GRO file. Must be a single line.
    atoms
        An instance of Structure containing the atoms to write.
    box
        The periodic box as a 3x3 matrix.
    """
    # Print the title
    print('TITLE ' + title, file=outfile)

    # Print the box
    print(pdbBoxString(box), file=outfile)

    # Print the atoms
    for idx, atname, resname, resid, x, y, z in atoms:
        print(pdbline % (idx % 1e5, atname, resname, "", resid % 1e5, '',
                         10*x, 10*y, 10*z, 0, 0, ''),
              file=outfile)


def write_summary(protein, membrane, solvent):
    charge  = protein.charge + membrane.charge
    plen = len(protein)
    print("; NDX Solute %d %d" % (1, protein and plen or 0), file=sys.stderr)
    print("; Charge of protein: %f" % protein.charge, file=sys.stderr)

    mlen = len(membrane)
    print("; NDX Membrane %d %d" % (1 + plen, membrane and plen + mlen or 0),
          file=sys.stderr)
    print("; Charge of membrane: %f" % membrane.charge, file=sys.stderr)
    print("; Total charge: %f" % charge, file=sys.stderr)

    slen = len(solvent)
    print("; NDX Solvent %d %d" % (1+plen+mlen, solvent and plen+mlen+slen or 0), file=sys.stderr)
    print("; NDX System %d %d" % (1, plen+mlen+slen), file=sys.stderr)
    print("; \"I mean, the good stuff is just INSANE\" --Julia Ormond",
          file=sys.stderr)


def write_top(outpath, molecules, title):
    """
    Write a basic TOP file.

    The topology is written in *outpath*. If *outpath* is en empty string, or
    anything for which ``bool(outpath) == False``, the topology is written on
    the standard error, and the header is omitted, and only what has been buit
    by Insane id displayed (e.g. Proteins are excluded).

    Parameters
    ----------
    outpath
        The path to the file to write. If empty, a simplify topology is
        written on stderr.
    molecules
        List of molecules with the number of them.
    title
        Title of the system.
    """
    topmolecules = []
    for i in molecules:
        if i[0].endswith('.o'):
            topmolecules.append(tuple([i[0][:-2]]+list(i[1:])))
        else:
            topmolecules.append(i)

    if outpath:
        # Write a rudimentary topology file
        with open(outpath, "w") as top:
            print('#include "martini.itp"\n', file=top)
            print('[ system ]', file=top)
            print('; name', file=top)
            print(title, file=top)
            print('\n', file=top)
            print('[ molecules ]', file=top)
            print('; name  number', file=top)
            print("\n".join("%-10s %7d"%i for i in topmolecules), file=top)
    else:
        # Here we only include molecules that have beed added by insane.
        # This is usually concatenated at the end of an existint top file.
        # As the existing file usually contain the proteins already, we do not
        # include them here.
        added_molecules = (molecule for molecule in topmolecules
                           if molecule[0] != 'Protein')
        print("\n".join("%-10s %7d"%i for i in added_molecules), file=sys.stderr)


def system_title(membrane, protein, lipids):
    (lipL, absL, relL), (lipU, absU, relU) = lipids
    if membrane.atoms:
        title  = "INSANE! Membrane UpperLeaflet>"+":".join(lipU)+"="+":".join([str(i) for i in relU])
        title += " LowerLeaflet>"+":".join(lipL)+"="+":".join([str(i) for i in relL])

        if protein:
            title = "Protein in " + title
    else:
        title = "Insanely solvated protein."
    return title


def write_structure(output, title, atoms, box):
    oStream = output and open(output, "w") or sys.stdout
    with oStream:
        if output.endswith(".gro"):
            write_gro(oStream, title, atoms, box.tolist())
        else:
            write_pdb(oStream, title, atoms, box.tolist())


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
