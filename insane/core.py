#!/usr/bin/env python3

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

import os
import sys
import random
import collections
import numpy as np

from simopt import opt_func

from . import lipids
from .pbc import PBC
from .structure import *
from .converters import *
from .constants import d2r
from ._data import SOLVENTS, CHARGES, APOLARS
from .options import OPTIONS


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


def _point(y, phi):
    r = np.sqrt(1-y*y)
    return np.cos(phi)*r, y, np.sin(phi)*r


def pointsOnSphere(n):
    return np.array([_point((2.*k+1)/n-1, k*2.3999632297286531) for k in range(n)])


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

    The PBC is changed **in place**.
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
            # We do not know what size the box should be. 
            # Let X and Y be the same.
            #T This does not agree with the default pbc being hexagonal...
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

    # Exclusion range for solvent
    excl = int(nz*options["solexcl"]/pbc.z)

    # Layer index at middle of box
    hz = int(0.5*nz) 

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
        for x, y, z in coord + options["lipradius"] * pointsOnSphere(options["lipdensity"]):
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
    solnames, solabs, solnums = list(zip(*solv))
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
        solnames.append("NA")
        num_sol.append(nna)
        solv.append("NA")
    if ncl:
        solnames.append("CL")
        num_sol.append(ncl)
        solv.append("CL")


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

    return sol, list(zip(solnames, num_sol))


def setup_membrane(pbc, protein, lipid, options):
    membrane = Structure()
    molecules = []

    lower, upper = lipid
    lipL, absL, relL = lower
    lipU, absU, relU = upper
    nabsL, nabsU = sum(absL), sum(absU)

    if not any((absL, relL, absU, relU)):
        return membrane, molecules

    # Update lipids name - add force field name to all lipids without 
    lipL = [lip if '.' in lip else options["forcefield"]+'.'+lip for lip in lipL]
    lipU = [lip if '.' in lip else options["forcefield"]+'.'+lip for lip in lipU]

    lo_lipd = np.sqrt(options["area"])
    up_lipd = np.sqrt(options["uparea"])

    #print('lipd', lo_lipd, up_lipd)

    num_up, num_lo = [], []

    # Lipids are added on grid positions, using the prototypes defined above.
    # If a grid position is already occupied by protein, the position is untagged.
    
    lipd = lo_lipd

    # Number of lipids in x and y in lower leaflet if there were no solute
    lo_lipids_x = int(pbc.x/lo_lipd+0.5)
    lo_lipids_y = int(pbc.y/lo_lipd+0.5)
    q = False
    while lo_lipids_x*lo_lipids_y < nabsL:
        if q:
            lo_lipids_x += 1
        else:
            lo_lipids_y += 1
        q = not q

    lo_lipdx    = pbc.x/lo_lipids_x
    lo_rlipx    = list(range(lo_lipids_x))
    lo_lipdy    = pbc.y/lo_lipids_y
    lo_rlipy    = list(range(lo_lipids_y))
    
    if options["uparea"]:
        lipd = up_lipd

    # Number of lipids in x and y in upper leaflet if there were no solute
    up_lipids_x = int(pbc.x/up_lipd+0.5)
    up_lipids_y = int(pbc.y/up_lipd+0.5)
    while up_lipids_x*up_lipids_y < nabsU:
        if q:
            up_lipids_x += 1
        else:
            up_lipids_y += 1
        q = not q

    up_lipdx    = pbc.x/up_lipids_x
    up_rlipx    = list(range(up_lipids_x))
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
    #print('->', len(grid_lo), len(grid_lo[0]))
    #print('->', len(grid_up), len(grid_up[0]))


    ## OLD GRID
        
    maxd = 1

    # If there is a protein, mark the corresponding cells as occupied
    if protein:
        # Extract the parts of the protein that are in either leaflet
        sphere = options["protradius"] * pointsOnSphere(options["protdensity"])

        # Calculate number density per cell
        mem_mask_lo = (0 > protein.coord[:,2]) & (protein.coord[:,2] > -2.4)
        prot_lo = protein.coord[mem_mask_lo, :]
        for i in prot_lo:
            for x, y, z in i + sphere:
                grid_lo[ int(lo_lipids_x*x/pbc.rx)%lo_lipids_x ][ int(lo_lipids_y*y/pbc.ry)%lo_lipids_y ] += 1

        mem_mask_up = (0 < protein.coord[:,2]) & (protein.coord[:,2] < 2.4)
        prot_up = protein.coord[mem_mask_up, :]
        for i in prot_up:
            for x, y, z in i + sphere:
                grid_up[ int(up_lipids_x*x/pbc.rx)%up_lipids_x ][ int(up_lipids_y*y/pbc.ry)%up_lipids_y ] += 1

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
    # Last, add lipids from command line, note first update the names with ff tag if needed 
    usrnames = [usrname if '.' in usrname else options["forcefield"]+'.'+usrname for usrname in options["lipnames"]]
    liplist.add_from_def(usrnames, options["lipheads"], options["liplinks"],
                         options["liptails"], options["lipcharge"])
    # Add lipid charges to global CHARGE index 
    for key in liplist:
        if liplist[key].charge != "0":
            CHARGES[key] = int(liplist[key].charge)

    if protein:
        resi = protein.atoms[-1][2]
    else:
        resi = 0

    inshift = options["indist"] / 2
    for leaflet, leaf_lip, lipd, lipdx, lipdy in [leaf_up, leaf_lo]:
        for lipid, pos in leaf_lip:
            # Increase the residue number by one
            resi += 1

            # Fetch the atom list with x, y, z coordinates
            try:
                at, ax, ay, az = zip(*liplist[lipid].build(diam=lipd))
            except KeyError as e:
                print(f"ERROR lipid name {e.args[0]} not found in database check lipids.dat, included mol files or specified definition strings")
                raise e 
            # The z-coordinates are spaced at 0.3 nm,
            # starting with the first bead at 0.15 nm
            az = [ leaflet*(inshift + (i-min(az)))*options["beaddist"] for i in az ]
            xx = np.array((ax, ay)).T

            # Set the random rotation for this lipid
            rangle   = 2*random.random()*np.pi
            if options["norotate"]:
                nx = pos[0] + lipdx/2 + [ kick*random.random() for i in az ]
                ny = pos[1] + lipdy/2 + [ kick*random.random() for i in az ]
            else:
                rcos     = np.cos(rangle)
                rsin     = np.sin(rangle)
                rcosx    = rcos*lipdx*2/3
                rcosy    = rcos*lipdy*2/3
                rsinx    = rsin*lipdx*2/3
                rsiny    = rsin*lipdy*2/3
                nx = np.dot(xx,(rcosx, -rsiny)) + pos[0] + lipdx/2 + [ kick*random.random() for i in az ]
                ny = np.dot(xx,(rsinx, rcosy)) + pos[1] + lipdy/2 + [ kick*random.random() for i in az ]

            # Add the atoms to the list
            memcoords.extend([(nx[i], ny[i], az[i]) for i in range(len(at))])
            mematoms.extend([(at[i], lipid, resi, 0, 0, 0) for i in range(len(at))])

    ##< Done building lipids

    membrane.coord = memcoords
    membrane.atoms = mematoms

    return membrane, molecules


@opt_func(OPTIONS)
def old_main(**options):

    molecules = []

    #########################################
    ## I. PROTEIN and other macromolecules ##
    #########################################

    # Read in the structures (if any)
    tm = [ Structure(i, options) for i in options["solute"] ]

    if tm:
        molecules.append(('Protein', len(tm)))

    xshifts  = [0] # Shift in x direction per protein


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

    #print(pbc.box)
    resize_pbc_for_lipids(pbc=pbc, relL=relL, relU=relU, absL=absL, absU=absU,
                          uparea=options["uparea"], area=options["area"],
                          hole=options["hole"], proteins=tm)
    #print(pbc.box)

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
        if '.' in i[0]:
            # Remove any -ff tags from molecules - WARNING no name can contain . as used as separator
            topmolecules.append(tuple([i[0].split('.')[1]]+list(i[1:])))
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
