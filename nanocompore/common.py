# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from warnings import warn

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

def file_header_contains (fn, field_names, sep="\t"):
    with open (fn, "r") as fp:
        header = fp.readline().rstrip().split (sep)
    for f in field_names:
        if not f in header:
            return False
    return True

def numeric_cast_list (l):
    """
    Cast str values to integer or float from a list
    """
    l2 = []
    for i in l:
        if type(i)== str:
            try:
                i = int(i)
            except ValueError:
                try:
                    i = float(i)
                except ValueError:
                    pass
        l2.append(i)
    return l2

def counter_to_str (c):
    """Transform a counter dict to a tabulated str"""
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m
