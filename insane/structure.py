
import numpy as np

from .converters import *
from ._data import SOLVENTS, CHARGES, APOLARS


def occupancy(grid, points, spacing=0.01):
    """Return a vector with the occupancy of each grid point for 
    given array of points"""
    distances = ((grid[:,None,:] - points[None,:,:])**2).sum(axis=2)
    occupied = (distances < spacing).sum(axis=1)
    return occupied


def isPDBAtom(l):
    return l.startswith("ATOM") or l.startswith("HETATM")

def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    ## ===>   atom name,   res name,     res id, chain,       x,            y,             z
    return (a[12:16], a[17:20], int(a[22:26]), a[21], float(a[30:38])/10, float(a[38:46])/10, float(a[46:54])/10)


# Reformatting of lines in structure file
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1"
pdbline = "ATOM  %5i %-4s %4s%1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s  "


def angle(a, b):
    p = (a*b).sum()
    q = np.sqrt((a**2).sum()*(b**2).sum())
    return np.arccos(np.clip(p/q, -1.0, 1.0))*180/np.pi 


def pdbBoxString(box):
    box = np.array(box)

    # Box vectors
    u, v, w  = box

    # Box vector lengths
    nu, nv, nw = np.sqrt((box**2).sum(axis=1))

    # Box vector angles
    alpha = nv*nw == 0 and 90 or angle(v, w)
    beta  = nu*nw == 0 and 90 or angle(u, w)
    gamma = nu*nv == 0 and 90 or angle(u, v)

    # Conversion from nm to A
    return pdbBoxLine % (10*nu, 10*nv, 10*nw, alpha, beta, gamma)


def groAtom(a):
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===>   atom name,   res name,     res id, chain,       x,          y,          z
    return (a[10:15], a[5:10],   int(a[:5]), " ", float(a[20:28]), float(a[28:36]), float(a[36:44]))


def groBoxRead(a):
    b = [float(i) for i in a.split()] + 6*[0] # Padding for rectangular boxes
    return b[0], b[3], b[4], b[5], b[1], b[6], b[7], b[8], b[2]


class Structure(object):
    def __init__(self, filename=None, options=None):
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

        if options:
            self.setup(**options)


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
            if '.' in resname:
                # Remove any -ff tags from molecules - WARNING no name can contain . as used as separator
                resname = resname.split('.')[1]
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
        if not mask.sum():
            # No cross section with membrane
            return 0
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
        m = self.coord.min(axis=0)
        r = self.coord.max(axis=0)-self.coord.min(axis=0)+1e-8

        # Number of grid cells
        n = (r/d + 0.5).astype('int')
        nx, ny, nz = n

        # Initialize grids
        atom = np.zeros(n+2)
        phobic = np.zeros(n+2)
        binned = (n * (self.coord - m) / r).astype('int')
        notdummy = np.array([ i[1] != "DUM" for i in self.atoms ])
        apolar = np.array([ i[1].strip() in APOLARS for i in self.atoms ])
        for i,j,k in binned[notdummy]:
            atom[i,j,k] += 1
        for i,j,k in binned[apolar]:
            phobic[i,j,k] += 1
        #atom = np.histogramdd(binned[notdummy], n+2)[0]

        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = (2 * phobic / atom) ** pw

        # Determine average density
        occupd = atom.astype('bool').sum()
        avdens = float(atom.sum())/occupd
        threshold = 0.1*avdens
        above = atom > threshold
        surface  = []
        #cgofile  = open('density.cgo', "w")
        #cgofile.write('[\n')
        for i,j,k in zip(*np.where(atom > threshold)):
            # Check the neighbouring cells;
            # if one of them is not occupied, count cell as surface
            if not (atom[i-1,j,k] and atom[i+1,j,k] and
                    atom[i,j-1,k] and atom[i,j+1,k] and
                    atom[i,j,k-1] and atom[i,j,k+1]):
                sx, sy, sz = m + (r*(i,j,k)+0.5*r)/n
                sw = ratio[i,j,k]
                surface.append((sx, sy, sz, sw))
                #cgofile.write("    7.0, %f, %f, %f, %f, \n"%(10*sx, 10*sy, 10*sz, 0.25*sw))
        #cgofile.write(']\n')
        #cgofile.close()

        surface = np.array(surface)
        sx, sy, sz, w = zip(*surface)
        W             = 1.0/sum(w)

        # Weighted center of apolar region; has to go to (0, 0, 0)
        sxm, sym, szm = np.average(surface[:,:-1], axis=0, weights=surface[:,-1])
        apolar_center = np.average(surface[:,:-1], axis=0, weights=surface[:,-1])

        # Place apolar center at origin
        self.center = -apolar_center

        # Determine weighted deviations from centers
        dev = (surface[:,:-1] - apolar_center) * surface[:, -1, None]

        # Covariance matrix for surface
        val, vec = np.linalg.eig(np.dot(dev.T, dev))
        vec = vec[:,val.argsort()[::-1]]

        # Rotate the coordinates
        self.coord = np.dot(self.coord, vec)

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

    def setup(self, **kwargs):
        # Center the protein and store the shift
        shift = self.center
        self.center = (0, 0, 0)

        ## 1. Orient with respect to membrane
        # Orient the protein according to the TM region, if requested
        # This doesn't actually work very well...
        if kwargs["orient"]:
            self.orient(kwargs["origriddist"], kwargs["oripower"])

        ## 4. Orient the protein in the xy-plane
        ## i. According to principal axes and unit cell
        self.rotate(kwargs["rotate"])

        ## 5. Determine the minimum and maximum x and y of the protein
        pmin, pmax = self.coord.min(axis=0), self.coord.max(axis=0)

        # At this point we should shift the subsequent proteins such
        # that they end up at the specified distance, in case we have
        # a number of them to do
        # y-shift is always -ycenter
        # x-shift is -xmin+distance+xmax(current)
        # xshifts.append(xshifts[-1]+pmax[0]+(options["distance"] or 0))

        ## 2. Shift of protein relative to the membrane center
        zshift = kwargs["memshift"]
        if not kwargs["center"]:
            zshift += shift[2]

        # Now we center the system in the rectangular
        # brick corresponding to the unit cell
        # If -center is given, also center z in plane
        self += (0, 0, zshift)

        # The z position is now correct with respect to the membrane
        # at z = 0. The x/y need to be set still


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
        print(pdbline % (idx % 1e5, atname[:4], resname[:3], "",
                         resid % 1e4, '', 10*x, 10*y, 10*z, 0, 0, ''),
              file=outfile)


def write_structure(output, title, atoms, box):
    oStream = output and open(output, "w") or sys.stdout
    # The title is always truncated to 80 characters to avoid 
    # Gromacs (or other programs) choking on them.
    with oStream:
        if output.endswith(".gro"):
            write_gro(oStream, title[:80], atoms, box.tolist())
        else:
            write_pdb(oStream, title[:80], atoms, box.tolist())

