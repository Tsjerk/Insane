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
import difflib
import functools
import os
import shutil
import shlex
import subprocess
import tempfile

from StringIO import StringIO

HERE = os.path.abspath(os.path.dirname(__file__))
INSANE = os.path.abspath(os.path.join(HERE, '../insane.py'))
DATA_DIR = os.path.join(HERE, 'data')
INSANE_SEED = '42'

# The arguments to test insane with are listed here. The tuple is used both to
# generate the references, and to run the tests.
# To add a test case, add the arguments to test in the tuple.
SIMPLE_TEST_CASES = (
    '-o test.gro',
    '-o test.gro -box 10,15,20',
    '-o test.gro -box 10,15,20 -l POPC',
    '-o test.gro -box 10,15,20 -sol W',
    '-o test.gro -box 10,15,20 -sol WF',
    '-o test.gro -box 10,15,20 -sol W -l POPC',
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
            o_index = i
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


def run_insane(arguments):
    command = [INSANE] + arguments
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               env={'INSANE_SEED': INSANE_SEED})
    out, err = process.communicate()
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


def run_and_compare(arguments, ref_gro, ref_stdout, ref_stderr):
    # Create the command as a list for subprocess.Popen.
    # The arguments can be pass to the current function as a string or as a
    # list of arguments. If they are passed as a string, they need to be
    # converted to a list.
    arguments =  _arguments_as_list(arguments)

    # The name of the output gro file must be provided to insane for insane to
    # work. Since we also need that file name, let's get it from insane's
    # arguments.
    gro_output = _gro_output_from_arguments(arguments)

    # We want insane to run in a temporary directory. This allows to keep the
    # file system clean, and it avoids mixing output of different tests.
    with tempdir():
        out, err = run_insane(arguments)
        assert os.path.exists(gro_output)
        compare(gro_output, ref_gro)
        compare(ContextStringIO(out), ref_stdout)
        compare(ContextStringIO(err), ref_stderr)


def test_simple_cases():
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    for case in SIMPLE_TEST_CASES:
        ref_gro = os.path.join(simple_case_ref_data, case + '.gro')
        ref_stdout = os.path.join(simple_case_ref_data, case + '.out')
        ref_stderr = os.path.join(simple_case_ref_data, case + '.err')
        # The test generator could yield run and compare directly. Bt, then,
        # the verbose display of nosetests gets crowded with the very long
        # names of the reference file, that are very redundant. Using a partial
        # function allows to have only the arguments for insane displayed.
        _test_case = functools.partial(
            run_and_compare,
            ref_gro=ref_gro,
            ref_stdout=ref_stdout,
            ref_stderr=ref_stderr)
        yield (_test_case, case)


def generate_simple_case_references():
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    for case in SIMPLE_TEST_CASES:
        arguments = _arguments_as_list(case)
        out_gro = _gro_output_from_arguments(arguments)
        ref_gro = os.path.join(simple_case_ref_data, case + '.gro')
        ref_stdout = os.path.join(simple_case_ref_data, case + '.out')
        ref_stderr = os.path.join(simple_case_ref_data, case + '.err')
        with tempdir():
            print(INSANE + ' ' + ' '.join(arguments))
            out, err = run_insane(arguments)
            with open(ref_stdout, 'w') as outfile:
                for line in out:
                    print(line, file=outfile, end='')
            with open(ref_stderr, 'w') as outfile:
                for line in err:
                    print(line, file=outfile, end='')
            shutil.copy2(out_gro, ref_gro)


if __name__ == '__main__':
    generate_simple_case_references()
