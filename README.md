# Insane - A simple, versatile tool for building coarse-grained simulation systems

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
    -sol W \
    -salt 0.15 \
    -o system.gro -p topol.top
```

## Installation

Insane is available on [pypi](https://pypi.python.org/pypi/insane) and
can be installed using `pip`.

```bash
pip install insane
```

We recommend the use of python _virtual environments_. The documentation for
[`venv`](https://docs.python.org/3/library/venv.html) is quite helpful, but many other methods for
using virtual environments are available.

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
be used with gromacs. If the `-p` argument is omitted, insane writes the
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
    -sol W
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
[martinize2](https://github.com/marrink-lab/vermouth-martinize) to prepare a
coarse grain protein structure.

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
    -l POPC -sol W \
    -f protein.gro -center
```

Here, the protein structure we give with the `-f` argument is centered in the
box along the Z axis with the `-center` argument.

With the `-pbc` option we set the shape of the periodic box to a prism with
a hexagonal base.

### Solvent mixtures

In the same way multiple lipids can be provided to set up complex lipid
bilayers, mixtures of solvent beads can be made.

This is particularly useful when building a Martini 2 system, which requires
solvation using normal (W) and anti-freeze waters (WF).

```bash
insane \
    -o system.gro -p topol.top \
    -ff M2 \
    -d 10 -l POPC \
    -f protein.gro -center \
    -sol W:90 -sol WF:10
```

This will set up a bilayer around a membrane protein and solvate the system with
90% normal waters, and with 10% anti-freeze waters.

Also note that the Martini 2 force field was selected using the `-ff M2` option.

### Selecting the lipid template force field

By default, insane will set up a Martini 3 system. This means that it draws from
Martini 3 lipid definitions. To change this, the `-ff M2` option can be used to
select Martini 2 lipids.

Lipids in lipids data files (e.g., `lipids.dat`) are prefixed by `M2.` or `M3.`
to indicate the force field version. To select a lipid with a particular force
field prefix directly, you can provide its fully specified name, which will
override the prefix given by the `-ff` option. For example, `M2.CHOL`,
`M3.CHOL`, and `M3d.CHOL` all refer to different cholesterol topologies. When
using `-ff M3` (the default value), you need to use the `M3d` prefix explicitly
to select that variant: `-l M3d.CHOL`.

### Adding custom lipid definitions

The default definitions of lipid templates are stored in the `insane/lipids.dat`
file.

A custom `lipids.dat` file can be added using the `-dat` option. For example, if
you have a custom lipids file called `custom.dat` that defines the special lipid
`ABCD`, you can use it as follows.

```bash
insane \
    -o system.gro \
    -dat custom.dat \
    -l ABCD -d 10
```

If the additional lipids file contains names that already exist, the previous definitions will be overwritten.

## Learn

Get usage information by running `insane -h`.

A number of tutorials using insane are available, including:

- Hands-on Martini 3 tutorial: [**Lipid bilayers: Complex mixtures**](https://cgmartini.nl/docs/tutorials/Martini3/LipidsII/)
- Book chapter: [**Building complex membranes with Martini 3**](https://doi.org/10.1016/bs.mie.2024.03.010)

  > Ozturk, T. N., König, M., Carpenter, T. S., Pedersen, K. B., Wassenaar, T. A., Ingólfsson, H. I.,
  > & Marrink, S. J. (2024). Building complex membranes with Martini 3. In _Methods in Enzymology_
  > (Vol. 701, pp. 237-285). Academic Press.
  > <https://doi.org/10.1016/bs.mie.2024.03.010>

## Contribute

Insane is hosted on [Github](https://github.com/Tsjerk/Insane). Please, report
any [issue](https://github.com/Tsjerk/Insane/issues) you encounter.

[gromacs]: http://www.gromacs.org
