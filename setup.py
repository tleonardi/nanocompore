#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup
import nanocompore as package

# Collect info in a dictionnary for setup.py
setup(
    name =  package.__name__,
    version = package.__version__,
    description = package.__description__,
    url = "https://github.com/tleonardi/nanocompore",
    author = 'Tommaso Leonardi and Adrien Leger',
    author_email = 'tom@tleo.io',
    license = 'GPLv3',
    python_requires ='>=3.5',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3'],
    install_requires = [
        'numpy>=1.14.0',
        'scipy>=1.1.0',
        'tqdm>=4.23.4',
        "pyfaidx>=0.5.4.1",
        "matplotlib>=2.2.2",
        "seaborn>=0.9.0",
        "pandas>=0.23.3",
        "statsmodels>=0.9.0",
        "scikit-learn>=0.20",
        "bedparse>=0.1.2",
        "pyyaml>=5.0"],
    packages = [package.__name__],
    package_data = {__name__: ["models/kmers_model_RNA_r9.4_180mv.tsv"]},
    entry_points = {'console_scripts': ['nanocompore=nanocompore.__main__:main']})
