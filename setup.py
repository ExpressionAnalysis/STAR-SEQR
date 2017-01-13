#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
from io import open
from sys import stderr


def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


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
    install_requires=['six', 'pandas>=0.18.0', 'pysam>=0.9.0', 'primer3-py', 'intervaltree_bio'],
    ext_modules=[libssw_ext],
    package_data={"starseqr_utils": ["resources/*"]},
    scripts=['starseqr.py', 'scripts/starseqr2bedpe.py'],
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'],
    keywords=['rna', 'rna-seq', 'fusions', 'chimeric', 'star']
)

stderr.write("Installation was successful!\n")
