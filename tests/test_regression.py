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
Regression tests for insane.

This test suite runs the insane command line with various set of arguments, and
assess that the results correspond to the result obtained with previous
versions.

Notice that these tests do not assess that the results are correct. Instead,
they assess that changes do not affect the behaviour of the program.

If ran as a script, this generate the reference files expected by the tests. If
ran usinf pytest or nosetest, this executes insane with a series of arguments
and compares the output to the reference.
"""

import contextlib
import functools
import glob
import mock
import os
import random
import shutil
import shlex
import subprocess
import sys
import tempfile
import textwrap

from nose.tools import assert_equal, assert_raises

import insane.cli
import utils

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest

HERE = os.path.abspath(os.path.dirname(__file__))
# INSANE = os.path.abspath(os.path.join(HERE, '../insane.py'))
INSANE = utils.which('insane')
DATA_DIR = os.path.join(HERE, 'data')
INPUT_DIR = os.path.join(HERE, 'data', 'inputs')
INSANE_SEED = '42'

# The arguments to test insane with are listed here. The tuple is used both to
# generate the references, and to run the tests.
# To add a test case, add the arguments to test in the tuple.
SIMPLE_TEST_CASES = [
    '-o test.gro',
    '-o test.gro -box 10,15,20',
    '-o test.gro -box 10,15,20 -l POPC',
    '-o test.gro -box 10,15,20 -sol W',
    '-o test.gro -box 10,15,20 -sol WF',
    '-o test.gro -box 10,15,20 -sol W -l POPC',
    '-o test.gro -l POPC -l DPPC -d 10',
    '-o test.gro -l POPC:2 -l DPPC:1 -d 10',
    ('-o test.gro -f CG1a0s.pdb', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -ring', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -orient', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -orient -od 0.2',
     '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -rotate princ', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -rotate 30', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -rotate 40', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -dm 3', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -center', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 20,30,40', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 20,30,40 -d 3',
     '1a0s'),
    ('-o test.pdb -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 20,30,40 -d 3',
     '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 2,3,4 -d 10',
     '1a0s'),
    ('-o test.gro -box 10,15,20 -l LOLO -alname LOLO '
     '-alhead C.P -allink G.G -altail CC.CDC'),
    ('-o test.gro -box 10,15,20 -l LOLO -l LOL2 '
     '-alname LOLO -alhead C.P -allink G.G -altail CC.CDC '
     '-alname LOL2 -alhead E.P -allink A.A -altail TCC.CDDC',
     None, 'multi-custom-lipids'),
    '-o test.pdb -box 10,15,20',
    ('-o test.pdb -f CG1a0s.pdb -p CG1a0s.top -l POPC -ring', '1a0s'),
    ('-o test.gro -pbc keep -f CG1a0s-box.pdb', '1a0s'),
    '-o test.gro -box 25',
    '-o test.gro -box 25,15,10',
    '-o test.gro -box 25,15,10,0,0,5,0,5,5',
    '-o test.gro -box 25,20,15,90,60,60',
    '-o bilayer.pdb -x 15 -y 15 -z 10 -l DOPC:0.5 -l DOPS:0.5 -sol W',
]

# Add test cases for all PBC options.
for pbc in ('hexagonal', 'rectangular', 'square', 'cubic', 'optimal'):
    SIMPLE_TEST_CASES.extend([
        ('-o test.gro -pbc {} -f CG1a0s.pdb -p CG1a0s.top'.format(pbc),
         '1a0s'),
        '-o test.gro -pbc {} -d 10'.format(pbc),
        '-o test.gro -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -hole 4 -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -disc 8 -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -disc 8 -hole 4 -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -hole 4 -pbc {} -d 10'.format(pbc),
        '-o test.gro -disc 8 -pbc {} -d 10'.format(pbc),
        '-o test.gro -disc 8 -hole 4 -pbc {} -d 10'.format(pbc),
    ])


def _arguments_as_list(arguments):
    """
    Return the arguments as a list as expected by subprocess.Popen.

    The arguments can be provided as a string that will be spitted to a list.
    They can also be provided as a list, then the list will be returned
    untouched.
    """
    try:
        arguments_list = shlex.split(arguments)
    except ValueError:
        arguments_list = arguments
    return arguments_list


def _split_case(case):
    """
    Get the arguments and the input directory from a test case.
    """
    if len(case) == 3:
        case_args, input_dir, alias = case
        if input_dir is not None:
            input_dir = os.path.join(INPUT_DIR, input_dir)
    elif len(case) == 2:
        case_args, input_dir = case
        input_dir = os.path.join(INPUT_DIR, input_dir)
        alias = case_args
    else:
        case_args = case
        input_dir = None
        alias = case_args
    return case_args, input_dir, alias


def _reference_path(arguments, alias=None):
    """
    Get the path to the reference files for the simple test cases.
    """
    arg_list = _arguments_as_list(arguments)

    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    base_name = arguments if alias is None else alias
    return os.path.join(simple_case_ref_data, base_name)


def _run_internal(arguments):
    os.environ['INSANE_SEED'] = INSANE_SEED
    random.seed(INSANE_SEED)
    command = [INSANE] + arguments
    out = StringIO()
    err = StringIO()
    with mock.patch('random.random', return_value=0.1):
        with utils._redirect_out_and_err(out, err):
            returncode = insane.cli.main(command)
    out = out.getvalue()
    err = err.getvalue()
    return out, err, returncode


def run_insane(arguments, input_directory=None, runner=_run_internal):
    """
    Run insane with the given arguments

    Insane is run in a copy of `input_directory`.
    """
    # Copy the content of the input directory in the current directory if an
    # input directory is provided.
    if input_directory is not None:
        for path in glob.glob(os.path.join(input_directory, '*')):
            if os.path.isdir(path):
                shutil.copytree(path, '.')
            else:
                shutil.copy2(path, '.')
    out, err, returncode = runner(arguments)
    print("** Insane exited with return code {}.".format(returncode))
    if returncode:
        print(err)
    return out, err, returncode


def compare(output, reference):
    """
    Assert that two files are identical.
    """
    out_file = utils._open_if_needed(output)
    ref_file = utils._open_if_needed(reference)
    with out_file, ref_file:
        lines_zip = zip(out_file, ref_file)
        for out_line, ref_line in lines_zip:
            assert_equal(out_line, ref_line)
        extra_out = list(out_file)
        extra_ref = list(ref_file)
        assert_equal(extra_out, [])
        assert_equal(extra_ref, [])


def compare_directories(directory, ref_directory, ignore=()):
    extra_files = []
    missing_files = []
    for root, dirs, files in os.walk(directory):
        for single_file in [f for f in files if f not in ignore]:
            path = os.path.join(root, single_file)
            rel_path = os.path.relpath(path, directory)
            ref_path = os.path.join(ref_directory, rel_path)
            if os.path.exists(ref_path):
                if os.path.splitext(ref_path)[-1] == '.gro':
                    utils.assert_gro_equal(path, ref_path)
                compare(path, ref_path)
            else:
                extra_files.append(rel_path)
    for root, dirs, files in os.walk(ref_directory):
        for single_file in [f for f in files if f not in ignore]:
            path = os.path.join(root, single_file)
            rel_path = os.path.relpath(path, ref_directory)
            ref_path = os.path.join(directory, rel_path)
            if not os.path.exists(ref_path):
                missing_files.append(rel_path)
    assert not bool(extra_files), ("The following files are unexpected: {}"
                                   .format(extra_files))
    assert not bool(missing_files), ("The following files are missing: {}"
                                     .format(missing_files))


def run_and_compare(arguments, input_dir, ref_dir, runner):
    """
    Run insane and compare its output against a reference
    """
    # Create the command as a list for subprocess.Popen.
    # The arguments can be pass to the current function as a string or as a
    # list of arguments. If they are passed as a string, they need to be
    # converted to a list.
    arguments = _arguments_as_list(arguments)

    ref_stdout = os.path.join(ref_dir, 'stdout')
    ref_stderr = os.path.join(ref_dir, 'stderr')

    # We want insane to run in a temporary directory. This allows to keep the
    # file system clean, and it avoids mixing output of different tests.
    with utils.tempdir():
        out, err, returncode = run_insane(arguments, input_dir, runner=runner)
        assert not returncode
        compare(utils.ContextStringIO(out), ref_stdout)
        compare(utils.ContextStringIO(err), ref_stderr)
        compare_directories('./', ref_dir,
                            ignore=('stderr', 'stdout', 'testlog'))


def test_simple_cases_internal():
    """
    This function generates test functions for nosetests. These test functions
    calls insane's main function with the argument listed in SIMPLE_TEST_CASES.
    """
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir, alias = _split_case(case)
        ref_dir = _reference_path(case_args, alias)
        # The test generator could yield run and compare directly. Bt, then,
        # the verbose display of nosetests gets crowded with the very long
        # names of the reference file, that are very redundant. Using a partial
        # function allows to have only the arguments for insane displayed.
        _test_case = functools.partial(
            run_and_compare,
            ref_dir=ref_dir,
            runner=_run_internal)
        yield (_test_case, case_args, input_dir)

class TestGroTester(object):
    """
    Test if the comparison of GRO file catches the differences.
    """
    ref_gro_content = """\
    INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
    4
        1POPC   NC3    1   2.111  14.647  11.951
        1POPC   PO4    2   2.177  14.644  11.651
        1POPC   GL1    3   2.128  14.642  11.351
        1POPC   GL2    4   1.961  14.651  11.351
    10 10 10"""

    def test_equal(self):
        """
        Make sure that identical files do not fail.
        """
        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            utils.assert_gro_equal('ref.gro', 'ref.gro')

    def test_diff_x(self):
        """
        Make sure that error in coordinates is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.353  # Is not within tolerance
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_in_tolerance(self):
        """
        Make sure that small errors in coordinates are not caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.352  # Is within tolerance
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            utils.assert_gro_equal('content.gro', 'ref.gro')

    def test_diff_natoms(self):
        """
        Make sure that differences in number of atom is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        6
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
            1POPC   C1A    5   2.125  14.651  11.051
            1POPC   D2A    6   2.134  14.602  10.751
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_title(self):
        """
        Make sure that a different title is caught.
        """
        gro_content = """\
        A different title
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_box(self):
        """
        Make sure that a different box is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1POPC   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
        10 9.9 10 9.08 4 54"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')

    def test_diff_field(self):
        """
        Make sure that a difference in a field is caught.
        """
        gro_content = """\
        INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1
        4
            1POPC   NC3    1   2.111  14.647  11.951
            1DIFF   PO4    2   2.177  14.644  11.651
            1POPC   GL1    3   2.128  14.642  11.351
            1POPC   GL2    4   1.961  14.651  11.351
        10 10 10"""

        with utils.tempdir():
            with open('ref.gro', 'w') as outfile:
                print(textwrap.dedent(self.ref_gro_content),
                      file=outfile, end='')
            with open('content.gro', 'w') as outfile:
                print(textwrap.dedent(gro_content), file=outfile, end='')
            assert_raises(AssertionError, utils.assert_gro_equal,
                          'content.gro', 'ref.gro')


def generate_simple_case_references():
    """
    Run insane to generate reference files for the simple regression tests.

    Run insane with the arguments listed in SIMPLE_TEST_CASES. The output
    files, the standard output, and the standard error are stored in the
    DATA_DIR/simple_case directory.
    """
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir, alias = _split_case(case)
        arguments = _arguments_as_list(case_args)
        ref_dir = _reference_path(case_args, alias)
        ref_stdout = os.path.join(ref_dir, 'stdout')
        ref_stderr = os.path.join(ref_dir, 'stderr')

        if os.path.exists(ref_dir):
            shutil.rmtree(ref_dir)
        os.mkdir(ref_dir)

        with utils.in_directory(ref_dir):
            print(INSANE + ' ' + ' '.join(arguments))
            out, err, _ = run_insane(arguments, input_dir)
            with open(ref_stdout, 'w') as outfile:
                for line in out:
                    print(line, file=outfile, end='')
            with open(ref_stderr, 'w') as outfile:
                for line in err:
                    print(line, file=outfile, end='')


def clean_simple_case_references():
    """
    Delete reference files for the simple tests if they are not in use anymore.
    """
    simple_test_cases = [_split_case(case)[2] for case in SIMPLE_TEST_CASES]
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    for path in glob.glob(os.path.join(simple_case_ref_data, '*')):
        base_name = os.path.basename(os.path.splitext(path)[0])
        if base_name not in simple_test_cases:
            print(path)
            os.remove(path)


def main():
    """
    Command line entry point.
    """
    help_ = """
Generate or clean the reference files for insane's regression tests.

{0} gen: generate the files
{0} clean: clean the unused files

nosetests -v: run the tests
""".format(sys.argv[0])
    commands = {'gen': generate_simple_case_references,
                'clean': clean_simple_case_references}
    if len(sys.argv) != 2:
        print(help_, file=sys.stderr)
        sys.exit(1)
    try:
        commands[sys.argv[1]]()
    except KeyError:
        print("Unrecognized keyword '{}'.".format(sys.argv[1]))
        print(help_, file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
