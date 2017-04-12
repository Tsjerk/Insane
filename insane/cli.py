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

import sys

import simopt
from . import core
from .options import OPTIONS


def main(argv):
    ## TEMPORARY ---
    # Exception to be defined in insane
    class InsaneBuildException(BaseException): pass
    ## <---

    ## OPTIONS
    # Parse options
    try:
        options = OPTIONS.parse(argv[1:])
    except simopt.SimoptHelp:
        print(OPTIONS.help(argv[1:]))
        return 0
    except simopt.MissingMandatoryError as e:
        print(e)
        return 3
    except simopt.Usage as e:
        print(e)
        return 1

    ## WORK
    try:
        system = core.insane(**options)
    except InsaneBuildException as e:
        print(e)
        return 2

    ## OUTPUT
    # Build atom list
    # Build topology
    # Build index



    (molecules,
     protein,
     membrane,
     solvent,
     lipids,
     box) = core.old_main(**options)

    title = core.system_title(membrane, protein, lipids)
    atoms = protein + membrane + solvent

    core.write_summary(protein, membrane, solvent)
    core.write_structure(
        output=options['output'],
        title=title,
        atoms=atoms,
        box=box,
    )
    core.write_top(options['topology'], molecules, title)

    return 0


def cli():
    sys.exit(main(sys.argv))
