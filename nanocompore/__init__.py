# -*- coding: utf-8 -*-
# https://github.com/python-poetry/poetry/pull/2366#issuecomment-652418094
try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)

__description__ = 'Software package that identifies raw signal changes between two conditions from https://github.com/jts/nanopolish resquiggled dRNA-Seq data.'
