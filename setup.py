#!/usr/bin/env python3

from __future__ import print_function, absolute_import
from setuptools import setup

# Read the version from a file to be sure it is consistent with the version
# in the package.
with open('VERSION.txt') as infile:
    version = infile.readline().strip()

setup(
    name='insane',
    version=version,

    description="A versatile tool for building membranes and/or solvent with proteins.",

    url='https://github.com/Tsjerk/Insane',

    # Author details
    author='Tsjerk A. Wassenaar',

    license='GPLv2',

    classifiers=[
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',

        'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',

        'Programming Language :: Python :: 3'
    ],

    install_requires=['numpy', 'simopt>=0.4.0'],

    tests_requires=['nose'],

    packages=['insane'],
    package_data={'insane': ['VERSION.txt', '*.dat', 'data/*.dat']},

    entry_points={
        'console_scripts': [
            'insane = insane.cli:cli',
        ],
    },

)
