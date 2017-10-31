
import numpy as np

RTOL = 1e-8


class PBCException(Exception):
    pass


class PBC(object):
    def __init__(self, shape=None, box=None, xyz=None, distance=None, membrane=None, protein=None, disc=None, hole=None):
        self.box = None

        # Validate the shape.
        if shape not in ('cubic', 'rectangular', 'square',
                         'hexagonal', 'optimal', 'keep', None):
            raise PBCException('"{}" is not a known PBC shape.'.format(shape))

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


