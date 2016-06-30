#!/usr/bin/env python

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

from __future__ import print_function

from nose.tools import assert_equal, assert_raises

import contextlib
import functools
import glob
import os
import shutil
import shlex
import subprocess
import tempfile

from StringIO import StringIO

HERE = os.path.abspath(os.path.dirname(__file__))
INSANE = os.path.abspath(os.path.join(HERE, '../insane.py'))
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
    '-o test.gro -pbc hexagonal -x 10 -y 15 -z 20',
    ('-o test.gro -f CG1a0s.pdb', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -ring', '1a0s'),
]

# Add test cases for all PBC options.
for pbc in ('hexagonal', 'rectangular', 'square', 'cubic', 'optimal'):
    SIMPLE_TEST_CASES.append(
        ('-o test.gro -pbc {} -f CG1a0s.pdb -p CG1a0s.top'.format(pbc), '1a0s')
    )


class ContextStringIO(StringIO):
    """
    StringIO but usable as a context manager.

    StringIO can not completely pass as a file, because it is not usable as a
    context manager in a 'with' statement. This class adds the context manager
    hability to StringIO. It does nothing when the context manager is either
    entered or exited, but it can be used in a 'with' statement.
    """
    def __enter__(self):
        """
        Does nothing when entering a 'with' statement.
        """
        pass

    def __exit__(self, *args):
        """
        Does nothing when exiting a 'with' statement.
        """
        pass


@contextlib.contextmanager
def tempdir():
    """
    Context manager that moves in a temporary directory.

    .. code::

        # We are in the initial working directory
        with tempdir():
            # We are now in a temporary directory
            ...
        # We are back to the initial working directory.
        # The temporary does not exist anymore.
    """
    return_dir = os.getcwd()
    dirpath = tempfile.mkdtemp()
    try:
        os.chdir(dirpath)
        yield dirpath
    finally:
        os.chdir(return_dir)
        shutil.rmtree(dirpath)


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


def _gro_output_from_arguments(arguments):
    """
    Find the file name of the GRO output provided as argument to insane.

    The file name is passed to insane via the '-o' argument. If the argument is
    provided several times, then only the last one is considered.

    This function reads the arguments provided as a list of arguments.
    """
    for i, argument in reversed(list(enumerate(arguments))):
        if argument == '-o':
            break
    else:
        raise ValueError('Output GRO name is not provided to insane. '
                         'Missing -o argument.')
    return arguments[i + 1]


def _open_if_needed(handle):
    """
    Return handle if it is a ContextStringIO instance else try to open it.
    """
    if isinstance(handle, ContextStringIO):
        return handle
    return open(handle)


def _split_case(case):
    """
    Get the arguments and the input directory from a test case.
    """
    if len(case) == 2:
        case_args, input_dir = case
        input_dir = os.path.join(INPUT_DIR, input_dir)
    else:
        case_args = case
        input_dir = None
    return case_args, input_dir


def _reference_path(arguments):
    """
    Get the path to the reference files for the simple test cases.
    """
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    ref_gro = os.path.join(simple_case_ref_data, arguments + '.gro')
    ref_stdout = os.path.join(simple_case_ref_data, arguments + '.out')
    ref_stderr = os.path.join(simple_case_ref_data, arguments + '.err')
    return ref_gro, ref_stdout, ref_stderr


def run_insane(arguments, input_directory=None):
    # Copy the content of the input directory in the current directory if an
    # input directory is provided.
    if input_directory is not None:
        for path in glob.glob(os.path.join(input_directory, '*')):
            if os.path.isdir(path):
                shutil.copytree(path, '.')
            else:
                shutil.copy2(path, '.')
    command = [INSANE] + arguments
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               env={'INSANE_SEED': INSANE_SEED})
    out, err = process.communicate()
    if process.returncode:
        print(err)
    return out, err


def compare(output, reference):
    """
    Assert that two files are identical.
    """
    out_file = _open_if_needed(output)
    ref_file = _open_if_needed(reference)
    with out_file, ref_file:
        for out_line, ref_line in zip(out_file, ref_file):
            assert_equal(out_line, ref_line)
        assert_raises(StopIteration, next, out_file)
        assert_raises(StopIteration, next, ref_file)


def run_and_compare(arguments, input_dir, ref_gro, ref_stdout, ref_stderr):
    # Create the command as a list for subprocess.Popen.
    # The arguments can be pass to the current function as a string or as a
    # list of arguments. If they are passed as a string, they need to be
    # converted to a list.
    arguments = _arguments_as_list(arguments)

    # The name of the output gro file must be provided to insane for insane to
    # work. Since we also need that file name, let's get it from insane's
    # arguments.
    gro_output = _gro_output_from_arguments(arguments)

    # We want insane to run in a temporary directory. This allows to keep the
    # file system clean, and it avoids mixing output of different tests.
    with tempdir():
        out, err = run_insane(arguments, input_dir)
        assert os.path.exists(gro_output)
        compare(gro_output, ref_gro)
        compare(ContextStringIO(out), ref_stdout)
        compare(ContextStringIO(err), ref_stderr)


# This function generates test functions for nosetests. These test functions
# execute insane with the argument listed in SIMPLE_TEST_CASES.
# Do not add a docstring to this function. The docstring would be displayed in
# the verbose mode of nosetests preventing to distinguish among the different
# tests that are generated.
def test_simple_cases():
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir = _split_case(case)
        ref_gro, ref_stdout, ref_stderr = _reference_path(case_args)
        # The test generator could yield run and compare directly. Bt, then,
        # the verbose display of nosetests gets crowded with the very long
        # names of the reference file, that are very redundant. Using a partial
        # function allows to have only the arguments for insane displayed.
        _test_case = functools.partial(
            run_and_compare,
            ref_gro=ref_gro,
            ref_stdout=ref_stdout,
            ref_stderr=ref_stderr)
        yield (_test_case, case_args, input_dir)


def generate_simple_case_references():
    """
    Run insane to generate reference files for the simple regression tests.

    Run insane with the arguments listed in SIMPLE_TEST_CASES. The output GRO
    file, the standard output, and the standard error are stored in the
    DATA_DIR/simple_case directory.
    """
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir = _split_case(case)
        arguments = _arguments_as_list(case_args)
        out_gro = _gro_output_from_arguments(arguments)
        ref_gro, ref_stdout, ref_stderr = _reference_path(case_args)
        with tempdir():
            print(INSANE + ' ' + ' '.join(arguments))
            out, err = run_insane(arguments, input_dir)
            with open(ref_stdout, 'w') as outfile:
                for line in out:
                    print(line, file=outfile, end='')
            with open(ref_stderr, 'w') as outfile:
                for line in err:
                    print(line, file=outfile, end='')
            shutil.copy2(out_gro, ref_gro)


if __name__ == '__main__':
    generate_simple_case_references()
