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
Test the functions from insane.utils.
"""

from insane import utils

def test_open_resource_root():
    # Test that a ressource at the root of the module can be opened
    assert list(utils.iter_resource('lipids.dat'))


def test_open_resource_nested():
    # Test that a ressource in a subdirectory of the module can be opened
    assert list(utils.iter_resource('data/diacylester.dat'))
