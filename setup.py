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
    version='0.0.4',
    description='',
    long_description=read('README.md'),
    license='Custom',
    author='Jeff Jasper',
    author_email='jasper1918@gmail.com',
    url='https://github.com/ExpressionAnalysis/STAR-SEQR',
    packages=find_packages(),
    install_requires=['pandas>=0.18.0', 'pysam==0.9.0', 'primer3-py', 'intervaltree_bio'],
    package_data={"starseqr_utils": ["resources/*"]},
    scripts=['starseqr.py'],
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'],
)

stderr.write("Installation successful!\n")

sys.exit(0)
