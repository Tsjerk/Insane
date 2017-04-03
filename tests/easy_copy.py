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
Nose plugin that rewrite the test description for better regression tests

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
    """
    Color insane arguments in regression tests
    """
    categories = (
        (('l', 'alname', 'alhead', 'allink', 'altail'), '1;33'),  # lipids
        (('box', 'pbc', 'x', 'y', 'z'), '38;5;214'),  # PDB
        (('f', ), '31'),  # solute
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
    """
    Rewrites insane regression test description for easy copy-paste
    """
    result = arguments
    input_dir = None
    func_names = ('test_regression.test_simple_cases',
                  'test_regression.test_simple_cases_internal')
    func_name = None
    for name in func_names:
        if arguments.startswith(name + '('):
            func_name = name
            break
    if func_name is not None:
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
    """
    Make regression test descriptions easier to read and use.
    """
    enabled = True

    def __init__(self, *args, **kwargs):
        super(ImprovedDisplay, self).__init__(*args, **kwargs)
        self.color_verbose = self.enabled
        self.easy_copy = self.enabled
        self.replace_input = self.enabled

    def options(self, parser, env):
        """
        Registers the commandline option, defaulting to enabled.
        """
        super(ImprovedDisplay, self).options(parser, env)
        parser.add_option(
            "--no-color", action="store_false",
            default=True, dest="color_verbose",
            help="Disable color in regression tests."
        )
        parser.add_option(
            "--no-easy-copy", action="store_false",
            default=True, dest="easy_copy",
            help="Disable easy copy paste of regression tests."
        )
        parser.add_option(
            "--np-replace-input", action="store_false",
            default=True, dest="replace_input",
            help=("Do not replace input files by their absolute path "
                  "in regression test descriptions.")
        )

    def configure(self, options, conf):
        """
        Configure plugin. Plugin is enabled by default.
        """
        super(ImprovedDisplay, self).configure(options, conf)
        self.color_verbose = options.color_verbose
        self.easy_copy = options.easy_copy
        self.replace_input = options.replace_input
        try:
            self.enabled = options.color_verbose or options.easy_copy
        except AttributeError:
            self.enabled = False

    def describeTest(self, test):
        """
        Nose hook that rewrites test descriptions.
        """
        args = test.id()
        if self.easy_copy:
            args = replace_function(args, self.replace_input)
        if self.color_verbose:
            args = color_arguments(args)
        return args

    @staticmethod
    def report(stream):
        """
        Nose hook that prints at the end of the tests.
        """
        stream.write('Color is beautiful\n')
