# Installation

## Create a clean virtual environment

Ideally, before installation, create a clean **python3.5+** virtual environment to deploy the package. **Python 2 is not supported**. For example one can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):

```bash
mkvirtualenv nanocompore -p python3.6
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -n nanocompore python=3.6
```

## Dependencies

Nanocompore relies on a the following robustly maintained third party python libraries:

* numpy >= 1.14.0
* scipy >= 1.1.0
* tqdm >= 4.23.4
* pyfaidx >= 0.5.4.1
* matplotlib >= 2.2.2
* seaborn >= 0.9.0
* pandas >= 0.23.3
* statsmodels >= 0.9.0
* scikit-learn >= 0.20
* bedparse >= 0.1.2

The correct versions of packages are installed together with the software when using pip.

## Option 1: Direct installation with pip from github (recommended)

* To install the package with an ssh key

```bash
pip3 install git+ssh://git@github.com/tleonardi/nanocompore.git
```

* To install the package with https/ssh

```bash
pip3 install git+https://github.com/tleonardi/nanocompore.git
```

## Option 2: Clone the repository and install locally in develop mode

With this option, the package will be locally installed in “editable” or “develop” mode. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to modify the code and/or participate to the development of the package (see [contribution guidelines](https://github.com/tleonardi/nanocompore/blob/master/CONTRIBUTING.md)).

```bash
# Clone repo localy
git clone https://github.com/tleonardi/nanocompore.git

# Enter in repo directory
cd nanocompore

# Make setup.py executable
chmod u+x setup.py

# Install with pip3
pip3 install -e ./
```
