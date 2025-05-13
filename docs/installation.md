# Installation

## Create a clean virtual environment

Ideally, before installation, create a clean **python3.10+** virtual environment to deploy the package. The recommended python version is 3.12. **Python 2 is not supported**. For example one can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):

```bash
mkvirtualenv nanocompore -p python3.12
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -n nanocompore python=3.12
```

## Option 1: Direct installation with pip from PyPi or conda (recommended)

```bash
pip3 install nanocompore
```

```bash
conda install -c bioconda nanocompore
```

## Option 2: Clone the repository and install locally in develop mode with uv

With this option, the package will be locally installed in *editable* or *develop mode*. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to modify the code and/or participate to the development of the package (see [contribution guidelines](contributing.md)).


Make sure you have **uv** installed ([instructions](https://docs.astral.sh/uv/getting-started/installation/#installation-methods)).

```bash
# Clone repo localy
git clone https://github.com/tleonardi/nanocompore.git

# Enter in repo directory
cd nanocompore

# Install the dependencies with uv
uv sync
```

## Testing the installation (Optional)

If `nanocompore` is installed in develop mode (Option 2) it is recommended to run the unit tests and integration tests to verify the installation. The test framework `pytest` needs to be installed manually to run the test in the virtual environment where `NanoCompore` is also installed.


Then run `pytest` in the directory where `nanocompore` was downloaded.
```bash
uv run pytest tests
```

If all the tests are successful you should get a similar output:

```text
====================== test session starts ==========================
platform linux -- Python 3.12.8, pytest-8.3.5, pluggy-1.5.0
rootdir: /home/mzdravkov/nanocompore
configfile: pyproject.toml
plugins: cov-6.1.1, xdist-3.6.1, jaxtyping-0.2.38
collected 51 items

tests/test_common.py .                                         [  1%]
tests/test_comparisons.py ...........                          [ 23%]
tests/test_config.py ..................                        [ 58%]
tests/test_eventalign_collapse.py ...                          [ 64%]
tests/test_main.py .                                           [ 66%]
tests/test_postprocessing.py .......                           [ 80%]
tests/test_preprocessing.py .                                  [ 82%]
tests/test_run.py .......                                      [ 96%]
tests/test_uncalled4.py ..                                     [100%]

===================== 51 passed, 7 warnings in 18.98s =================
```

