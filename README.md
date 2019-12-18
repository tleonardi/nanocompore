![Nanocompore](docs/pictures/Nanocompore_logo.png)

---

[![GitHub license](https://img.shields.io/github/license/tleonardi/nanocompore)](https://github.com/tleonardi/nanocompore/blob/master/LICENSE)
[![Python](https://img.shields.io/badge/Python-%3E%3D3.6-yellow)](https://www.python.org/)
[![BioRxiv](https://img.shields.io/badge/BioRxiv-10.1101%2F843136%20-red)](https://www.biorxiv.org/content/10.1101/843136v1.full)

[![PyPI version](https://badge.fury.io/py/nanocompore.svg)](https://badge.fury.io/py/nanocompore)
[![Downloads](https://pepy.tech/badge/nanocompore)](https://pepy.tech/project/nanocompore)
[![Build Status](https://travis-ci.com/tleonardi/nanocompore.svg?token=2uTrW9fP9RypfMALjksc&branch=master)](https://travis-ci.com/tleonardi/nanocompore)

---

**Nanocompore identifies differences in ONT nanopore sequencing raw signal corresponding to RNA modifications by comparing 2 samples**

Nanocompore compares 2 ONT nanopore direct RNA sequencing datasets from different experimental conditions expected to have a significant impact on RNA modifications. It is recommended to have at least 2 replicates per condition. For example one can use a control condition with a significantly reduced number of modifications such as a cell line for which a modification writing enzyme was knocked-down or knocked-out. Alternatively, on a smaller scale transcripts of interests could be synthesized in-vitro.

**Full documentation is available at http://nanocompore.rna.rocks**

[![Nanocompore](docs/pictures/worflow.png)](http://nanocompore.rna.rocks)

## Companion repositories

* [NanoCompore_pipeline](https://github.com/tleonardi/nanocompore_pipeline): Nextflow pipeline to preprocess data for NanoCompore
* [Nanocompore_analysis](https://github.com/tleonardi/nanocompore_paper_analyses): Analyses performed with Nanocompore for the BioRxiv preprint
* [NanopolishComp](https://github.com/tleonardi/NanopolishComp): Collapse Nanopolish eventalign output per kmer, required before running NanoCompore

## Main authors

* Tommaso Leonardi - tom {at} tleo.io
* Adrien Leger - aleg {at} ebi.ac.uk
