# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from nanocompore.NanocomporeError import NanocomporeError
from tqdm import tqdm, tqdm_notebook

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#

def mkdir (fn):
    """Create directory recursivelly. Raise IO error if path exist or if error at creation
    """
    if os.path.isdir (fn):
        raise NanocomporeError ("The output folder specified already exists")
    else:
        try:
            os.makedirs (fn)
        except:
            raise NanocomporeError ("Error creating output folder {}".format(fn))

def access_file (fn, **kwargs):
    """Check if the file is readable
    """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)

def counter_to_str (c):
    """Transform a counter dict to a tabulated str"""
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m

def mytqdm (**kwargs):
    try:
        if get_ipython().__class__.__name__ == 'ZMQInteractiveShell':
            return tqdm_notebook(**kwargs)
        else:
            return tqdm(**kwargs)
    except NameError:
        return tqdm(**kwargs)
