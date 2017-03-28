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
import sys
import tempfile

from StringIO import StringIO

def which(program):
    """Determine full path of executable *program* on :envvar:`PATH`.

    (Jay at http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python)
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        real_program = realpath(program)
        if is_exe(real_program):
            return real_program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

HERE = os.path.abspath(os.path.dirname(__file__))
#INSANE = os.path.abspath(os.path.join(HERE, '../insane.py'))
#INSANE = '/home/jon/Envs/tsjerk/bin/insane'
INSANE = which('insane')
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
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -orient -od 0.2', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -rotate princ', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -rotate 30', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -rotate 40', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -dm 3', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -center', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 20,30,40', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 20,30,40 -d 3', '1a0s'),
    ('-o test.gro -f CG1a0s.pdb -p CG1a0s.top -l POPC -box 2,3,4 -d 10', '1a0s'),
    '-o test.gro -box 10,15,20 -l LOLO -alname LOLO -alhead C.P -allink G.G -altail CC.CDC',
    ('-o test.gro -box 10,15,20 -l LOLO -l LOL2 '
     '-alname LOLO -alhead C.P -allink G.G -altail CC.CDC '
     '-alname LOL2 -alhead E.P -allink A.A -altail TCC.CDDC', None, 'multi-custom-lipids'),
    '-o test.pdb -box 10,15,20',
    ('-o test.pdb -f CG1a0s.pdb -p CG1a0s.top -l POPC -ring', '1a0s'),
    ('-o test.gro -pbc keep -f CG1a0s-box.pdb', '1a0s'),
    '-o test.gro -box 25',
    '-o test.gro -box 25,15,10',
    '-o test.gro -box 25,15,10,0,0,5,0,5,5',
    '-o test.gro -box 25,20,15,90,60,60',
]

# Add test cases for all PBC options.
for pbc in ('hexagonal', 'rectangular', 'square', 'cubic', 'optimal'):
    SIMPLE_TEST_CASES.extend([
        ('-o test.gro -pbc {} -f CG1a0s.pdb -p CG1a0s.top'.format(pbc), '1a0s'),
        '-o test.gro -pbc {} -d 10'.format(pbc),
        '-o test.gro -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -hole 4 -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -disc 8 -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -disc 8 -hole 4 -pbc {} -d 10 -dz 5'.format(pbc),
        '-o test.gro -hole 4 -pbc {} -d 10'.format(pbc),
        '-o test.gro -disc 8 -pbc {} -d 10'.format(pbc),
        '-o test.gro -disc 8 -hole 4 -pbc {} -d 10'.format(pbc),
    ])


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
    out_struct = _gro_output_from_arguments(_arguments_as_list(arguments))
    out_format = os.path.splitext(out_struct)[-1]
    simple_case_ref_data = os.path.join(DATA_DIR, 'simple_case')
    base_name = arguments if alias is None else alias
    ref_gro = os.path.join(simple_case_ref_data, base_name + out_format)
    ref_stdout = os.path.join(simple_case_ref_data, base_name + '.out')
    ref_stderr = os.path.join(simple_case_ref_data, base_name + '.err')
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
        extra_out = list(out_file)
        extra_ref = list(ref_file)
        assert_equal(extra_out, [])
        assert_equal(extra_ref, [])


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
        case_args, input_dir, alias = _split_case(case)
        ref_gro, ref_stdout, ref_stderr = _reference_path(case_args, alias)
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
    for case in SIMPLE_TEST_CASES:
        case_args, input_dir, alias = _split_case(case)
        arguments = _arguments_as_list(case_args)
        out_gro = _gro_output_from_arguments(arguments)
        ref_gro, ref_stdout, ref_stderr = _reference_path(case_args, alias)
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
    help_ = """
Generate or clean the reference files for insane's regression tests.

{0} gen: generate the files
{0} clean: clean the unused files

nosetests -v: run the tests
""".format(sys.argv[0])
    commands = {'gen': generate_simple_case_references,
                'clean': clean_simple_case_references,}
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
