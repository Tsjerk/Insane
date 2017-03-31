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
Nose plugin that rewrite the test description for easier to use regression tests

The plugin allows to:

* color the arguments of the test cases to identify what they refer to;
* write the test case descriptionis as insane commands that can be copy-pasted;
* when writing the test case description as insane commands, replace the input
  file names by absolute paths.
"""

import os
import re
import shlex
from nose.plugins.base import Plugin


def color_arguments(arguments):
    categories = (
        (('l', 'alname', 'alhead', 'allink', 'altail'), '1;33'),  # lipids
        (('box', 'pbc', 'x', 'y', 'z'), '38;5;214'),  # PDB
        (('f', ),  '31'),  # solute
        (('d', 'dz'), '34'),  # distance to the box
        (('sol', ), '36'),  # solvent
    )
    template_from = "(-({}) [a-zA-Z0-9:.,/-]+)(\\s|$|')"
    template_to = '\\033[{}m\\1\\033[0m\\3'
    for option, color in categories:
        arguments = re.sub(template_from.format('|'.join(option)),
                           template_to.format(color),
                           arguments)
    return arguments 


def replace_function(arguments, replace_files=True):
    result = arguments
    input_dir = None
    func_name = 'test_regression.test_simple_cases'
    if arguments.startswith(func_name):
        reduced_name = arguments[len(func_name) + 1:-1]
        insane_args, input_dir = shlex.split(reduced_name)
        if replace_files:
            arg_list = shlex.split(insane_args)
            new_list = []
            while arg_list:
                arg = arg_list.pop(0)
                new_list.append(arg)
                if arg == '-f':
                    arg = arg_list.pop(0)
                    new_list.append(os.path.join(input_dir, arg))
            insane_args = ' '.join(new_list)
        result = 'insane ' + insane_args[:-1]
    return result


class ImprovedDisplay(Plugin):
    enabled = True
    name = 'color-verbose'
    score = 1600

    def options(self, parser, env):
        """Registers the commandline option, defaulting to enabled.
        """
        parser.add_option(
            "--no-color".format(self.name), action="store_false",
            default=True, dest="color_verbose",
            help="Disable color."
        )
        parser.add_option(
            "--no-easy-copy".format(self.name), action="store_false",
            default=True, dest="easy_copy",
            help="Disable easy copy paste."
        )
        parser.add_option(
            "--np-replace-input".format(self.name), action="store_false",
            default=True, dest="replace_input",
            help="Do not replace input files by their absolute path."
        )

    def configure(self, options, conf):
        """Configure plugin. Plugin is enabled by default.
        """
        self.color_verbose = options.color_verbose
        self.easy_copy = options.easy_copy
        self.replace_input = options.replace_input
        self.config = conf
        try:
            self.enabled = options.color_verbose or options.easy_copy
        except AttributeError:
            self.enabled = False

    def describeTest(self, test):
        args = test.id()
        if self.easy_copy:
            args = replace_function(args, self.replace_input)
        if self.color_verbose:
            args = color_arguments(args)
        return args

    def report(self, stream):
        stream.write('Color is beautiful\n')
