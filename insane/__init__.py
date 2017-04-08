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

__authors__ = ["Tsjerk A. Wassenaar", "Jonathan Barnoud"]
__year__    = "2017"


# Read the version from a file to make sure 
# that it is consistent with the one in setup.py.
import os
here = os.path.dirname(__file__)
try:
    with open(os.path.join(here, 'VERSION.txt')) as infile:
        __version__ = infile.readline().strip()
except:
    __version__ = 'unknown'
# Avoid poluting the namespace with useless stuff.
del os
del here


from .core import *
