# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from warnings import warn
from collections import *

#~~~~~~~~~~~~~~CUSTOM EXCEPTION CLASS~~~~~~~~~~~~~~#
class NanocomporeError (Exception):
    """ Basic exception class for nanocompore module """
    pass

class NanocomporeWarning (Warning):
    """ Basic Warning class for nanocompore module """
    pass

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

def numeric_cast_list (l):
    """
    Cast str values to integer or float from a list
    """
    l2 = []
    for v in l:
        l2.append(numeric_cast(v))
    return l2

def numeric_cast_dict (keys, values):
    """
    Cast str values to integer or float from a list
    """
    d = OrderedDict()
    for k, v in zip(keys, values):
        d[k] = numeric_cast(v)
    return d

def numeric_cast (v):
    if type(v)== str:
        try:
            v = int(v)
        except ValueError:
            try:
                v = float(v)
            except ValueError:
                pass
    return v

def counter_to_str (c):
    """Transform a counter dict to a tabulated str"""
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m

def all_values_in (required_val_list, all_val_list):
    """"""
    for v in required_val_list:
        if not v in all_val_list:
            return False
    return True
