#Insane - A simple, versatile tool for building coarse-grained simulation systems

Insane (INSert membrANE) is a versatile tool to build coarse-grained simulation
systems containing solutes, lipid bilayers, and/or solvents. It is initially
aimed at the Martini force field but can be used for other models. It can
write PDB or gromacs files.

Using insane, you can build a system containing a transmembrane protein
included in a lipid bilayers, with water and ions with the following command
line:

```bash
insane \
    -f cg1a0s.pdb \
    -l POPC \
    -sol W:90 -sol WF:10 \
    -salt 0.8 \
    -o system.gro -p topol.top
```

## Installation

Insane can be installed either as a single-file program, or as a python module.

### Install insane as a single-file program

Insane depends on python 2.7 and on the numpy python library. Make sure you
have these two requirements installed:

```bash
python2.7 -c 'import numpy'
```

If this command line outputs an error message, you need to instal either python
or numpy.

Download the latest version of insane as a single executable.

Run insane: `./insane -h`.

### Install insane as a python module

Insane is available on pypi and can be installed using pip. To install insane
for the current user, run

```bash
pip2 install --user insane
```

Pip installs insane in `~/.local/lib/python2.7/site-packages`. Check that the
directory exists and that it is in your `PYTHON_PATH` and in your `PATH`.

To install insane system wide, run

```bash
pip2 install insane
```

We recommend the use of python virtual environments. Read more about them on
the MDAnalysis website.

## Quick start

### Create a simple lipid bilayer

The following command creates a hydrated POPC bilayer:

```bash
insane -l POPC -box 10,10,7 -sol W -o bilayer.gro -p topol.top
```

The `-l` argument tells insane to include the given lipid. The number of lipids
to include is determined from the area per lipid and the size of the box; the
default area per lipid in 0.6nm².

The `-box` argument defines the dimensions of the box. Here, we set an
orthorhombic periodic box that is 10 nm along the X and Y axes, and 7 nm along
the Z axis.

The `-sol` argument tells insane to include some solvent in the box. Here we
include Martini water. The total amount is defined from the free volume of the
box, and the solvent diameter (0.5 nm per default).

The `-o` and `-p` arguments specify the output of the program. With `-o
bilayer.gro`, insane will write the system structure in a GRO file usable with
gromacs. Insane can also write PDB files if the output file as the '.pdb' file
extension. With `-p topol.top`, insane will write a template TOP file to be
used with gromacs. If the `-p` argument is omitted, insane writes the content
of the system on the standard output in a form that can be appended at the end
of an existing TOP file.

### Create a more complex bilayer

Insane can build lipid mixtures. Here we build a ternary mixture containing
a fully saturated lipid (DPPC), a polyunsaturated lipid
(dilinoleyl-phosphatidylcholine, called DIPC in the Martini force field), and
cholesterol.

```bash
insane \
    -o bilayer.gro -p topol.top \
    -x 20 -x 20 -z 10 \
    -l DPPC:4 -l DIPC:3 -l CHOL:3 \
    -sol W:90 -sol W:10
```

The periodic box can be defined in multiple ways. Here we build an orthorhombic
box by setting the X, Y, and Z axes separately with `-x`, `-y`, and `-z`,
respectively.

The `-l` option can be provided multiple times to setup multiple lipid types.
Here we setup DPPC, DIPC, and cholesterol with a 4:3:3 ratio.

The `-sol` option can also be provided multiple times. Here the solvent is
composed of 90% of regular Martini water and 10% Martini anti-freeze water.

### Solvate a protein

Insane can setup protein systems. Here we create a box of solvent around
a protein. Note that insane does not prepare the protein itself. Use martinize
to coarse grain a protein structure.

```bash
insane -o system.gro -p topol.top -f protein.pdb -d 7 -sol W -salt 0 -exclude 0
```
Solutes, including proteins, can be read from a GRO or a PDB file provided with
the `-f` argument.

With the `-d 7` argument, insane build the periodic box such that the distance
between two periodic images is greater than 7 nm. By default, in the absence of
a membrane, insane build a rhombic dodecahedral box. Setting the periodic box
with `-box`, or setting the box geometry with `-pbc`, overwrite that default.


### Insert a protein in a membrane

### Create a box of solvent

## Get help

## Contribute

