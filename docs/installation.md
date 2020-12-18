# Installation

## Create a clean virtual environment

Ideally, before installation, create a clean **python3.6+** virtual environment to deploy the package. **Python 2 is not supported**. For example one can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):

```bash
mkvirtualenv nanocompore -p python3.7
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -n nanocompore python=3.7
```

## Dependencies

Nanocompore relies on a the following robustly maintained third party python libraries:

* python = ">=3.6.1"
* numpy = "~1.19"
* scipy = "~1.5"
* tqdm = "~4"
* pyfaidx = "~0.5"
* matplotlib = "~3.1"
* seaborn = "~0"
* pandas = "~0.25"
* statsmodels = "~0.12"
* scikit-learn = "~0.23"
* bedparse = "~0.2"
* pyyaml = "~5"
* loguru = "~0.5"

The correct versions of packages are installed together with the software when using pip.

## Option 1: Direct installation with pip from PyPi or conda (recommended)

```bash
pip3 install nanocompore
```

```bash
conda install -c bioconda nanocompore
```

## Option 2: Clone the repository and install locally in develop mode with poetry

With this option, the package will be locally installed in *editable* or *develop mode*. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to modify the code and/or participate to the development of the package (see [contribution guidelines](contributing.md)).

```bash
# Clone repo localy
git clone https://github.com/tleonardi/nanocompore.git

# Enter in repo directory
cd nanocompore

# Install with pip3
pip install poetry
poetry install
```

## Testing the installation (Optional)

If `nanocompore` is installed in develop mode (Option 2) it is recommended to run the unit tests and integration tests to verify the installation. The test framework `pytest` needs to be installed manually to run the test in the virtual environment where `NanoCompore` is also installed.

```bash
pip install pytest
```

Then run `pytest` in the directory where `nanocompore` was downloaded.
```bash
pytest
```

If all the tests are successful you should get a similar output:

```text
====================================== test session starts =======================================
platform linux -- Python 3.6.6, pytest-4.3.0, py-1.7.0, pluggy-0.8.1
rootdir: /home/aleg/Programming/nanocompore, inifile:
collected 40 items                                                                               

tests/test_Integration.py .................                                                [ 42%]
tests/test_SampCompDB.py ....                                                              [ 52%]
tests/test_TxComp.py ...............                                                       [ 90%]
tests/test_Whitelist.py ....                                                               [100%]

================================== 40 passed in 160.40 seconds ===================================
```
