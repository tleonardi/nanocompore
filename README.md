# Nanocompore

---

**Software package that identifies raw signal changes between two conditions from [https://github.com/jts/nanopolish](nanopolish) resquiggled dRNA-Seq data.**

---

# Installation

Ideally, before installation, create a clean python3 virtual environment to deploy the package, using virtualenvwrapper for example ([see tutorial](http://www.simononsoftware.com/virtualenv-tutorial-part-2/)).

## Dependencies

Nanocompore relies on a few Python third party libraries (tqdm and numpy). The correct versions of packages are installed together with the software when using pip.

## Installation with pip

* To install the package

    ```pip3 install git+https://github.com/a-slide/nanocompore.git```

* To update the package:

    ```pip3 install git+https://github.com/a-slide/nanocompore.git --upgrade```

# Usage

Nanocompore compares 2 nanopore sequencing datasets from 2 experimental conditions expected to hava an impact on DNA/RNA modifications. In particular the software was optimized to work with direct RNA sequencing data.

To use nanocompore raw nanopore sequencing data have to be prepared using a basecaller (Albacore or MInKNOW. Ask ONT), your favorite long read aligner such as [Minimap2](https://github.com/lh3/minimap2) and resquiggled using [Nanopolish](https://github.com/jts/nanopolish). Finally eventalign data have to be collapsed per kmer using [NanopolishComp Eventalign_collapse](https://github.com/a-slide/NanopolishComp)

The preparation of data and the package usage are detailed in the [usage jupyter notebook](https://github.com/a-slide/nanocompore/blob/master/tests/nanocompore_usage.ipynb)

# Note to power-users and developers

Please be aware that **nanocompore** is an experimental package that is still under development. It was tested under Linux Ubuntu 16.04 and in an HPC environment running under Red Hat Enterprise 7.1.

You are welcome to contribute by requesting additional functionalities, reporting bugs or by forking and submitting pull requests

Thank you

