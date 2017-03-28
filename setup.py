#!/usr/bin/env python

from __future__ import print_function, absolute_import
from setuptools import setup, find_packages

setup(
    name='insane',
    version='1.0-dev',

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

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],

    install_require=['numpy'],

    tests_require=['nose'],

    packages=['insane'],

    entry_points={
        'console_scripts': [
            'insane = insane.cli:cli',
        ],
    },

)
