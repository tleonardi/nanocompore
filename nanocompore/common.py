# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from collections import *


#~~~~~~~~~~~~~~CUSTOM EXCEPTION CLASS~~~~~~~~~~~~~~#
class NanocomporeError (Exception):
    """ Basic exception class for nanocompore module """
    pass

class NanocomporeWarning (Warning):
    """ Basic Warning class for nanocompore module """
    pass

#~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~#
def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise NanocomporeError ("Error creating output folder `{}`".format(fn))

def access_file (fn, **kwargs):
    """ Check if the file is readable """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)

def numeric_cast_list (l):
    """ Cast str values to integer or float from a list """
    l2 = []
    for v in l:
        l2.append(numeric_cast(v))
    return l2

def numeric_cast_dict (keys, values):
    """ Cast str values to integer or float from a list """
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
    """ Transform a counter dict to a tabulated str """
    m = ""
    for i, j in c.most_common():
        m += "\t{}: {:,}".format(i, j)
    return m

def all_values_in (required_val_list, all_val_list):
    """Check that all values in required_val_list are found in all_val_list"""
    for v in required_val_list:
        if not v in all_val_list:
            return False
    return True

def jhelp (f:"python function or method"):
    """
    Display a Markdown pretty help message for functions and class methods (default __init__ is a class is passed)
    jhelp also display default values and type annotations if available.
    Undocumented options are not displayed.
    The docstring synthax should follow the markdown formated convention below
    * f
        Function or method to display the help message for
    """
    # For some reason signature is not always importable. In these cases the build-in help is called instead
    try:
        from IPython.core.display import display, Markdown, HTML
        import inspect
    except (NameError, ImportError) as E:
        NanocomporeWarning ("jupyter notebook is required to use this function. Please verify your dependencies")
        help(f)
        return

    if inspect.isclass(f):
        f = f.__init__

    if inspect.isfunction(f) or inspect.ismethod(f):

        # Parse arguments default values and annotations
        sig_dict = OrderedDict()
        for name, p in inspect.signature(f).parameters.items():
            sig_dict[p.name] = []
            # Get Annotation
            if p.annotation != inspect._empty:
                sig_dict[p.name].append(": {}".format(p.annotation))
            # Get default value if available
            if p.default == inspect._empty:
                sig_dict[p.name].append("(required)")
            else:
                sig_dict[p.name].append("(default = {})".format(p.default))

        # Parse the docstring
        doc_dict = OrderedDict()
        descr = []
        lab=None
        for l in inspect.getdoc(f).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    doc_dict[lab] = []
                elif lab:
                    doc_dict[lab].append(l)
                else:
                    descr.append(l)

        # Reformat collected information in Markdown synthax
        s = "---\n\n**{}.{}**\n\n{}\n\n---\n\n".format(f.__module__, f.__name__, " ".join(descr))
        for k, v in doc_dict.items():
            s+="* **{}** *{}*\n\n{}\n\n".format(k, " ".join(sig_dict[k]), " ".join(v))

        # Display in Jupyter
        display (Markdown(s))
