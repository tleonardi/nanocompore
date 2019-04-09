# -*- coding: utf-8 -*-

# Define self package variable
__version__ = "1.0.0b7"
description = 'Software package that identifies raw signal changes between two conditions from https://github.com/jts/nanopolish resquiggled dRNA-Seq data.'

# Collect info in a dictionary for setup.py
setup_dict = {
    "name": __name__,
    "version": __version__,
    "description": description,
    "url": "https://github.com/tleonardi/nanocompore",
    "author": 'Tommaso Leonardi and Adrien Leger',
    "author_email": 'tom {at} tleo.io / aleg {at} ebi.ac.uk',
    "license": "MIT",
   python_requires":'>=3.5',
    "classifiers": [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT',
        'Programming Language :: Python :: 3'],
    "install_requires": [
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
    "packages": [__name__],
    "package_data":{__name__: ["models/kmers_model_RNA_r9.4_180mv.tsv"]},
    "entry_points":{'console_scripts': ['nanocompore=nanocompore.nanocompore_main:main']}}
