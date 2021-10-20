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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301, USA.

"""
Test PBC related functions.
"""

import insane


def test_valid_pbc_shape():
    def try_shape(shape):
        try:
            insane.pbc.PBC(shape=shape, xyz=(10, 10, 10))
        except insane.pbc.PBCException as e:
            if 'is not a known PBC shape' in str(e):
                raise AssertionError()
            raise e

    shapes = ('cubic', 'rectangular', 'square', 'hexagonal', 'optimal', 'keep', None)
    for shape in shapes:
        yield try_shape, shape
    

def test_invalid_pbc_shape():
    def try_shape(shape):
        try:
            insane.pbc.PBC(shape=shape, xyz=(10, 10, 10))
        except insane.pbc.PBCException as e:
            if not 'is not a known PBC shape' in str(e):
                raise AssertionError()
        else:
            raise AssertError()

        shapes = ('invalid', '', 9, True, False)
        for shape in shapes:
            yield try_shape, shape
            
