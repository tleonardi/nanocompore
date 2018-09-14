# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from nanocompore.NanocomporeError import NanocomporeError

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
    if not os.path.isfile (fn):
        raise NanocomporeError("{} not found)".format (fn))
    if not os.access (fn, os.R_OK):
        raise NanocomporeError("{} not readable)".format (fn))

def counter_to_str (c):
    """Transform a counter dict to a tabulated str"""
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m
