#!/usr/bin/env python

from setuptools import setup, find_packages
import sys
import os
from sys import stderr
from ez_setup import use_setuptools
use_setuptools()


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='STAR-SEQR',
    version='0.0.1',
    description='',
    long_description=read('README.md'),
    license='APACHE',
    author='Jeff Jasper',
    author_email='jasper1918@gmail.com',
    url='https://github.com/jasper1918/STAR-SEQR',
    packages=find_packages(),
    install_requires=['cython', 'pandas>=0.18.0', 'pysam==0.9.0', 'primer3-py', 'intervaltree_bio'],
    package_data={"starseqr_utils": ["resources/*"]},
    scripts=['starseqr.py'],
    zip_safe=False,
)

stderr.write("Installation successful!\n")

sys.exit(0)
