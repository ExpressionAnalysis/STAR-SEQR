#!/usr/bin/env python

from setuptools import setup, find_packages
import sys
from sys import stderr
from ez_setup import use_setuptools
use_setuptools()

setup(
    description='',
    author='Jeff Jasper',
    url='https://github.com/jasper1918/STAR-SEQR',
    author_email='jasper1918@gmail.com',
    version='0.0.1',
    install_requires=['pandas>=0.18.0', 'pysam>=0.9.0', 'primer3-py', 'intervaltree_bio'],
    packages=['starseqr'],
    # packages=find_packages(),
    package_data={"starseqr": ["resources/*"]},
    ext_modules=[],
    scripts=['scripts/starseqr.py'],
    name='starseqr',
    license='APACHE',
    include='README.md',
)

stderr.write("Installation successful!\n")

sys.exit(0)
