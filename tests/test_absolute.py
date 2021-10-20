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
Test that a user can ask for an absolute number of lipids.
"""

import sys
from nose.tools import assert_equal, assert_raises
import utils
import insane

def test_absolute_lipid_only():
    with utils._redirect_out_and_err(sys.stdout, sys.stdout):
        (molecules,
         protein,
         membrane,
         solvent,
         lipids,
         box) = insane.core.old_main(lower=[('DOPC', 72, 0)],
                                     output='plop.gro', zdistance=15)
    npo4 = len([atom for atom in membrane.atoms if atom[0] == 'PO4'])
    assert_equal(npo4, 72 * 2)
    assert_equal(molecules, [('DOPC', 72), ('DOPC', 72)])
    assert_equal(lipids, ((('DOPC',), (72,), (0,)), (('DOPC',), (72,), (0,))))
    assert_equal(len(protein), 0)
    assert_equal(len(solvent), 0)


def test_absolute_lipid_assym():
    with utils._redirect_out_and_err(sys.stdout, sys.stdout):
        (molecules,
         protein,
         membrane,
         solvent,
         lipids,
         box) = insane.core.old_main(lower=[('DOPC', 72, 0)],
                                     upper=[('POPC', 89, 0)],
                                     output='plop.gro', zdistance=15)
    npo4 = len([atom for atom in membrane.atoms if atom[0] == 'PO4'])
    assert_equal(npo4, 72 + 89)
    assert_equal(molecules, [('POPC', 89), ('DOPC', 72)])
    assert_equal(lipids, ((('DOPC',), (72,), (0,)), (('POPC',), (89,), (0,))))
    assert_equal(len(protein), 0)
    assert_equal(len(solvent), 0)


def test_absolute_lipid_multi():
    with utils._redirect_out_and_err(sys.stdout, sys.stdout):
        (molecules,
         protein,
         membrane,
         solvent,
         lipids,
         box) = insane.core.old_main(lower=[('DOPC', 72, 0), ('POPC', 89, 0)],
                                     output='plop.gro', zdistance=15)
    npo4 = len([atom for atom in membrane.atoms if atom[0] == 'PO4'])
    assert_equal(npo4, (72 + 89) * 2)
    assert_equal(molecules, [('DOPC', 72), ('POPC', 89),
                             ('DOPC', 72), ('POPC', 89)])
    assert_equal(lipids, ((('DOPC', 'POPC'), (72, 89), (0, 0)),
                          (('DOPC', 'POPC'), (72, 89), (0, 0))))
    assert_equal(len(protein), 0)
    assert_equal(len(solvent), 0)
