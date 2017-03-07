#!/usr/bin/env python

import os
from setuptools import setup, find_packages, Extension
from io import open
from sys import stderr


def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths


su_tests = package_files('starseqr_utils/tests')

libssw_ext = Extension('_libssw_ext', sources=['ssw/src/ssw.c'], include_dirs=['ssw/src/'])

setup(
    name='starseqr',
    version=get_version(open('starseqr_utils/__init__.py', encoding='utf-8').read()),
    description='RNA-Fusion Calling with STAR',
    long_description=open('README.rst', encoding='utf-8').read(),
    license='Custom',
    author='Jeff Jasper',
    author_email='jasper1918@gmail.com',
    url='https://github.com/ExpressionAnalysis/STAR-SEQR',
    packages=find_packages(),
    install_requires=['cython', 'six', 'networkx','pandas >= 0.18.0', 'pysam >= 0.9.0', 'primer3-py', 'intervaltree_bio'],
    ext_modules=[libssw_ext],
    package_data={"starseqr_utils": ["resources/*"], '': su_tests},
    scripts=['starseqr.py'],
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'],
    keywords=['rna', 'rna-seq', 'fusions', 'chimeric', 'star'],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Environment :: Console",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: Implementation :: CPython",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ]
)

stderr.write("Installation was successful!\n")
