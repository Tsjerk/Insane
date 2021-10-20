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

import copy
from collections import namedtuple
import contextlib
import math
import os
import shutil
import sys
import tempfile

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

try:
    from itertools import izip_longest as zip_longest
except ImportError:
    from itertools import zip_longest

# GRO file format description. The key is the name of the field, the value is a
# tuple from which the first element is the first (included) and last
# (excluded) indices of the field in the line, and the second element the type
# of the field content.
GRO_FIELDS = {
    "resid": ((0, 5), int),
    "resname": ((5, 10), str),
    "name": ((10, 15), str),
    "index": ((15, 20), int),
    "x": ((20, 28), float),
    "y": ((28, 36), float),
    "z": ((36, 44), float),
}
GRO_TEMPLATE = ('{resid:>5}{resname:<5}{name:>5}{index:>5}'
                '{x:8.3f}{y:8.3f}{z:8.3f}')

GroDiff = namedtuple('GroDiff', 'linenum line ref_line fields')


class ContextStringIO(StringIO):
    """
    StringIO but usable as a context manager.

    StringIO can not completely pass as a file, because it is not usable as a
    context manager in a 'with' statement. This class adds the context manager
    ability to StringIO. It does nothing when the context manager is either
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
def in_directory(dirpath):
    return_dir = os.getcwd()
    try:
        os.chdir(dirpath)
        yield dirpath
    finally:
        os.chdir(return_dir)


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
    dirpath = tempfile.mkdtemp()
    try:
        with in_directory(dirpath):
            yield dirpath
    finally:
        shutil.rmtree(dirpath)


# realpath and which are copied from MDAnalysis.
# MDAnalysis is released under the GPL v2 license.
# Read the full license at
# <https://github.com/MDAnalysis/mdanalysis/blob/develop/LICENSE>
def realpath(*args):
    """Join all args and return the real path, rooted at /.
    Expands '~', '~user', and environment variables such as :envvar`$HOME`.
    Returns ``None`` if any of the args is ``None``.
    """
    if None in args:
        return None
    return os.path.realpath(
        os.path.expanduser(os.path.expandvars(os.path.join(*args)))
    )


def which(program):
    """Determine full path of executable *program* on :envvar:`PATH`.

    (Jay at http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python)
    """

    def is_exe(fpath):
        """
        Returns True is the path points to an executable file.
        """
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
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


def read_gro(stream):
    """
    Parse a gro file

    Read an iterable over the lines of a GRO file. Returns the title, the
    atoms, and the box. The atoms are returned as a list of dict, where each
    dict represents an atom. The keys of the atom dicts are

    * 'resid' for the residue number;
    * 'resname' for the residue name;
    * 'index' for the atom index as written in the file;
    * 'name' for the atom name;
    * 'x', 'y', and 'z' for the atom coordinates.

    The box is returned as a list of floats.

    .. note::

       The function does not read velocities. Also, it does not support
       variable precision.
    """
    # The two first lines are the header. The first line is the title of the
    # structure, the second line is the number of atoms.
    title = next(stream)
    natoms = int(next(stream))

    # Read the atoms according to the format described in GRO_FIELDS.
    # We stop when we reached the number of atoms declared.
    atoms = []
    for atom_count, line in enumerate(stream, start=0):
        if atom_count == natoms:
            break
        atoms.append({key: convert(line[begin:end].strip())
                      for key, ((begin, end), convert) in GRO_FIELDS.items()})
    else:
        raise ValueError('Box missing or invalid number of atoms declared.')

    # Read the box.
    box = [float(x) for x in line.split()]

    # Make sure there is nothing after the box.
    try:
        next(stream)
    except StopIteration:
        pass
    else:
        raise ValueError('Extra lines after the box or '
                         'invalid number of atoms declared')

    return title, atoms, box


def compare_gro(stream, ref_stream, tolerance=0.001):
    """
    Compare two gro files with a tolerance on the coordinates

    The `stream` and `ref_stream` arguments are iterable over the lines of two
    GRO files to compare. The tolerance is provided in nanometer so that

        abs(x1 - x2) <= tolerance

    The function returns a list of differences. Each difference is represented
    as a 'GroDiff' named tupple with the following field:

    * 'linenum': the number of the line where the difference occurs, the line
      count starts at 1 to make the difference easier to finc in a text editor;
    * 'line' and 'ref_line': the line that differ as it is written in 'stream'
      and in 'ref_stream', respectivelly;
    * 'fields': the specific field that differ beween the lines.

    The lines in the differences are re-formated from the information that have
    been parsed. The exact string may differ, especially if velocities were
    provided.

    The fields in the difference list are named after the GRO_FIELDS dict. It
    is `None` if the difference occurs in the title or the number of atoms. It
    is 'box' if for the box.

    If the number of atoms differ between the two files, then the atoms that
    exist in only one of the files are all counted as different and the field
    is set to `None`. If the box also differ, then its line number in the
    repport will be wrong for one of the files.
    """
    differences = []
    title, atoms, box = read_gro(stream)
    title_ref, atoms_ref, box_ref = read_gro(ref_stream)

    # Compare the headers
    if title != title_ref:
        diff = GroDiff(linenum=1, line=title, ref_line=title_ref, fields=None)
        differences.append(diff)
    if len(atoms) != len(atoms_ref):
        diff = GroDiff(linenum=2,
                       line=str(len(atoms)),
                       ref_line=str(len(atoms_ref)),
                       fields=None)
        differences.append(diff)

    # Compare the atoms
    atom_iter = enumerate(
        zip_longest(atoms, atoms_ref, fillvalue={}),
        start=3
    )
    for linenum, (atom, atom_ref) in atom_iter:
        if atom and atom_ref:
            # Both the atom and the reference atoms are defined.
            # We compare the fields.
            diff_fields = []
            for gro_field, (_, gro_type) in GRO_FIELDS.items():
                if gro_type == float:
                    # The field is a float, we compare with a tolerance.
                    error = math.fabs(atom[gro_field] - atom_ref[gro_field])
                    if error > tolerance:
                        diff_fields.append(gro_field)
                else:
                    # The field is an int or a string, we check equality.
                    if atom[gro_field] != atom_ref[gro_field]:
                        diff_fields.append(gro_field)
        else:
            # At least one of the atoms is not defined. They are counted as
            # different anyway, no need to compare anything.
            diff_fields.append(None)
        if diff_fields:
            # We found a difference, add it to the list.
            line = ref_line = ''
            if atom:
                line = GRO_TEMPLATE.format(**atom)
            if atom_ref:
                ref_line = GRO_TEMPLATE.format(**atom_ref)
            diff = GroDiff(linenum=linenum,
                           line=line,
                           ref_line=ref_line,
                           fields=diff_fields)
            differences.append(diff)

    # Compare the box
    if box != box_ref:
        diff = GroDiff(linenum=linenum + 1,
                       line=' '.join(map(str, box)),
                       ref_line=' '.join(map(str, box_ref)),
                       fields='box')
        differences.append(diff)

    return differences


def format_gro_diff(differences, outstream=sys.stdout, max_atoms=10):
    """
    Format differences between GRO files in a human readable way.
    """
    if not differences:
        # Do not display anything if the two gro files are identical.
        return

    # We do not want to modify the input
    differences = copy.copy(differences)

    # Display the differences in metadata first
    if differences and differences[0].linenum == 0:
        diff = differences.pop(0)
        print('The title is different:', file=outstream)
        print(diff.line, file=outstream)
        print(diff.ref_line, file=outstream)
    if differences and differences[0].linenum == 1:
        diff = differences.pop(0)
        print('The number of atoms is different! '
              '"{}" instead of "{}".'.format(diff.line, diff.ref_line),
              file=outstream)
    if differences and differences[-1].fields == 'box':
        diff = differences.pop(-1)
        print('The box is different:', file=outstream)
        print(diff.line, file=outstream)
        print(diff.ref_line, file=outstream)

    # Then display the atoms. Only display 'max_atoms' ones.
    if len(differences) > max_atoms:
        print('There are {} atoms that differ. '
              'Only displaying the {} first ones.'
              .format(len(differences), max_atoms),
              file=outstream)
    for diff in differences[:max_atoms]:
        print('On line {}:'.format(diff.linenum), file=outstream)
        print(diff.line, file=outstream)
        print(diff.ref_line, file=outstream)


def assert_gro_equal(path, ref_path):
    """
    Raise an AssertionError if two GRO files are not semantically identical.
    """
    diff_out = StringIO()
    with open(path) as stream, open(ref_path) as ref_stream:
        differences = compare_gro(stream, ref_stream)
        format_gro_diff(differences, outstream=diff_out)
        assert len(differences) == 0, '\n' + diff_out.getvalue()

def _open_if_needed(handle):
    """
    Return handle if it is a ContextStringIO instance else try to open it.
    """
    if isinstance(handle, ContextStringIO):
        return handle
    return open(handle)


@contextlib.contextmanager
def _redirect_out_and_err(stdout, stderr):
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    try:
        sys.stdout = stdout
        sys.stderr = stderr
        yield
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
