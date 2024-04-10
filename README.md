# Insane - A simple, versatile tool for building coarse-grained simulation systems

[![Build Status](https://travis-ci.org/Tsjerk/Insane.svg?branch=master)](https://travis-ci.org/Tsjerk/Insane)

Insane (INSert membrANE) is a versatile tool to build coarse-grained simulation
systems containing solutes, lipid bilayers, and/or solvents. It is initially
aimed at the [Martini force field](http://cgmartini.nl) but can be used for
other models. It can write PDB or [gromacs][] files.

Using insane, you can build a system containing a transmembrane protein
included in a lipid bilayer, with water and ions with the following command
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

Insane depends on python (either version 2 or 3) and on the numpy python library. Make sure you
have these two requirements installed:

```bash
python -c 'import numpy'
```

If this command line outputs an error message, you need to instal either python
or numpy.

Download the [latest version of
insane](https://github.com/Tsjerk/Insane/releases/download/v1.0rc1/insane) as
a single executable.

Run insane: `./insane -h`.

### Install insane as a python module

Insane is available on [pypi](https://pypi.python.org/pypi/insane/1.0rc1) and
can be installed using pip. To install insane for the current user, run

```bash
pip install --user insane
```

Pip installs the insane in `~/.local/lib/python<version>/site-packages`,
where `<version>` is the version of python. Check that the
directory exists and that it is in your `PYTHONPATH`. The insane program is
installed in `~/.local/bin`, make sure this directory is in your `PATH`.

To install insane system wide, run

```bash
sudo pip install insane
```

We recommend the use of python virtual environments. Read more about them on
the [MDAnalysis website](http://www.mdanalysis.org/2017/04/07/environments/).

## Quick start

### Create a simple lipid bilayer

The following command creates a hydrated POPC bilayer:

```bash
insane -l POPC -d 10 -dz 7 -sol W -o bilayer.gro -p topol.top
```

The `-l` argument tells insane to include the given lipid. The number of lipids
to include is determined from the area per lipid and the size of the box; the
default area per lipid is 0.6 nm².

With `-d` and `-dz`, we specify the dimensions of the box. Here we create
a hexagonal prism box where images are separated by 10 nm in the XY plane, and
7 nm in the Z dimension. The shape of the box can be changed with the `-pbc`
argument.

The `-sol` argument tells insane to include some solvent in the box. Here we
include Martini water. The total amount is defined from the free volume of the
box, and the solvent diameter (0.5 nm per default).

The `-o` and `-p` arguments specify the output of the program. With `-o
bilayer.gro`, insane will write the system structure in a [GRO
file](http://manual.gromacs.org/current/online/gro.html) usable with
[gromacs][]. Insane can also write PDB files if the output file has the '.pdb'
file extension. With `-p topol.top`, insane will write a template TOP file to
be used with [gromacs][]. If the `-p` argument is omitted, insane writes the
content of the system on the standard output in a form that can be appended at
the end of an existing TOP file.

### Create a more complex bilayer

Insane can build lipid mixtures. Here we build a ternary mixture containing
a fully saturated lipid (DPPC), a polyunsaturated lipid
(dilinoleyl-phosphatidylcholine, called DIPC in the Martini force field), and
cholesterol.

```bash
insane \
    -o bilayer.gro -p topol.top \
    -x 20 -y 20 -z 10 \
    -l DPPC:4 -l DIPC:3 -l CHOL:3 \
    -sol W:90 -sol WF:10
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
a protein. Note that insane does not prepare the protein itself. Use
[martinize](http://cgmartini.nl/index.php/tools2/proteins-and-bilayers) to
coarse grain a protein structure.

```bash
insane -o system.gro -p topol.top -f protein.pdb -d 7 -sol W -salt 0
```
Solutes, including proteins, can be read from a GRO or a PDB file provided with
the `-f` argument.

With the `-d 7` argument, insane builds the periodic box such that the distance
between two periodic images is greater than 7 nm. By default, in the absence of
a membrane, insane builds a rhombic dodecahedral box. Setting the periodic box
with `-box`, or setting the box geometry with `-pbc`, overwrite that default.

If the `-salt` argument is set, insane will add ions to the system. With
`-salt` set to 0, insane will add enough chloride or sodium ions to
neutralize the system charge. If `-salt` is set to any positive value, it is
read as the concentration of ions in molar.

### Insert a protein in a membrane

Now that we have seen how to build a membrane and how to build a protein
system, we can insert a protein in a lipid bilayer.

```bash
insane \
    -o system.gro -p topol.top \
    -d 10 -pbc hexagonal \
    -l POPC -sol PW \
    -f protein.gro -center
```

Here, the protein structure we give with the `-f` argument is centered in the
box along the Z axis with the `-center` argument.

We build this system with [Martini polarizable
water](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000810)
using `-sol PW`.

With the `-pbc` option we set the shape of the periodic box to a prism with
a hexagonal base.

### Changing templates 

Templates for both Martini 2 and 3 molecules are predefined within insane, use the `-ff M2` or `-ff M3` to switch between them. Specific template versions can also be specified directly in the lipid name e.g. use `-l M3.POPC` instead of `-l POPC`. 

Definitions of lipid templates are stored in the insane/lipids.dat file and can be viewed and edited there. 

## Get help

Get the list of all the arguments by running `insane -h`.

You can get additional help in the ["Tools" section of the Martini
forum](http://cgmartini.nl/index.php/component/kunena/9-tools).

As well as access tutorials using insane in [Martini 2 membrane tutorials](http://www.cgmartini.nl/index.php/tutorials-general-introduction-gmx5/bilayers-gmx5) and [Martini 3 membrane tutorials](https://doi.org/10.1016/bs.mie.2024.03.010).

## Contribute

Insane is hosted on [Github](https://github.com/Tsjerk/Insane). Please, report
there any [issue](https://github.com/Tsjerk/Insane/issues) you encounter.

[gromacs]: http://www.gromacs.org
