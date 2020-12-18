# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from collections import *
import gzip
from glob import iglob
from loguru import logger

#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#

class SuperParser ():
    def __init__ (self,
        fn,
        select_colnames=[],
        n_lines=None,
        sep="\t",
        comment="#",
        change_colnames={},
        cast_colnames={}):
        """
        Open a parser for field delimited file and return an iterator yield lines as namedtuples
        Transparently parse gziped and multiple file with the same header
        * fn
            Path to a field delimited file or a list of files or a regex or a mix of everything (Can be gzipped)
        * select_colnames
            List of column names to use parse and return
        * sep
            field separator
        * comment
            skip any line starting with this string
        * cast_colnames
            Dict corresponding to fields (based on colnames) to cast in a given python type
        """

        # Init logger and counter
        self.counter = Counter()

        # Save self variables
        self._sep = sep
        self._comment = comment
        self._n_lines = n_lines

        # Input file opening
        self.f_list = self._open_files (fn)

        # Define colnames based on file header. It needs to be the same for all the files to parse
        fn, fp = self.f_list[0]
        logger.debug("Reading header from file: {}".format(fn))
        self.colnames = self._get_first_line_header(fp)
        self.select_idx = []

        # Find selected colnames if needed
        if select_colnames:
            for field in select_colnames:
                try:
                    idx = self.colnames.index(field)
                    self.select_idx.append(idx)
                except ValueError:
                    raise SuperParserError ("Required field `{}` not found in files `{}`".format(field, fn))
            self.colnames = select_colnames

        # Check header for subsequent files if any
        if len(self.f_list)>1:
            for fn, fp in self.f_list[1:]:
                colnames = self._get_first_line_header(fp)
                if self.select_idx:
                    for idx, col in zip(self.select_idx, self.colnames):
                        if idx >= len(colnames) or colnames[idx] != col:
                            raise SuperParserError ("Inconsistant headers between `{}` and `{}`".format(self.f_list[0][0], fn))
                else:
                    if colnames != self.colnames:
                        raise SuperParserError ("Inconsistant headers between `{}` and `{}`".format(self.f_list[0][0], fn))

        # Change colnames is required
        if change_colnames:
            for i, col in enumerate(self.colnames):
                if col in change_colnames:
                    self.colnames[i] = change_colnames[col]

        logger.debug("Column names from header: `{}`".format(" / ".join(self.colnames)))

        # Save initial number of columns, custom namedtuple and dtype index
        self.ncols = len(self.colnames)
        self.lt = namedtuple("lt", self.colnames)
        self.cast_colnames = cast_colnames

    #~~~~~~~~~~~~~~MAGIC AND PROPERTY METHODS~~~~~~~~~~~~~~#

    def __len__ (self):
        size = 0
        for fn, fp in self.f_list:
            size+= int(os.path.getsize(fn))
        return size-self._header_len

    def __enter__ (self):
        return self

    def close (self):
        for i, j in self.counter.most_common():
            logger.debug("{}: {}".format(i, j))
        for fn, fp in self.f_list:
            try:
                logger.debug ("Closing file:{}".format(fn))
                fp.close()
            except Exception as E:
                logger.warning (E)

    def __exit__(self, exception_type, exception_val, trace):
        self.close()

    def __iter__ (self):
        # Iterate over files
        for fn, fp in self.f_list:
            logger.debug("Starting to parse file {}".format(fn))
            # Iterate over line in file
            for line in fp:
                self.counter["Lines Parsed"]+=1

                if self._comment and line.startswith(self._comment):
                    self.counter["Comment lines skipped"]+=1
                    continue
                try:
                    line = self._parse_line(line)
                    self.counter["Lines successfully parsed"]+=1
                    yield line
                    # early stopping condition
                    if self._n_lines and self.counter["Lines successfully parsed"] == self._n_lines:
                        return
                except (SuperParserError, TypeError) as E:
                    self.counter["Malformed or Invalid Lines"]+=1
                    logger.debug(E)
                    logger.debug("File {}: Invalid line {}".format(fn, line))
            logger.debug("End of file: {}".format(fn))
        logger.debug("All files done")

    #~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~#

    def _get_first_line_header (self, fp):
        header_line = next(fp)
        header_list = header_line.rstrip().split(self._sep)
        return header_list

    def _parse_line (self, line):

        # Split line
        line = line.rstrip().split(self._sep)

        # Select field if needed
        if self.select_idx:
            line = [line[i] for i in self.select_idx]

        # Raise error if the length of the line is inconsistent with the header
        if len(line) != self.ncols:
            raise SuperParserError("Invalid Number of fields found")
        line_d = OrderedDict(zip(self.colnames,line))

        # Cast values according to provided types of lambda functions
        if self.cast_colnames:
            for colname, cast in self.cast_colnames.items():
                try:
                    line_d[colname] = cast(line_d[colname])
                except Exception:
                    raise SuperParserError("Cannot cast field in required type")

        # Return parsed line as a dict
        return line_d

    def _open_files (self, fn_list):
        """Transparently open files, lists, regex, gzipped or not"""
        f_list = []

        # Standard input
        if fn_list is 0:
            fn = "stdin"
            fp = open(0)
            return [(fn,fp)]

        # Cast single file or regex to list
        if isinstance(fn_list, str):
            fn_list = [fn_list]

        if isinstance(fn_list, (list, tuple, set)):
            for fn_regex in fn_list:
                for fn in iglob(fn_regex):
                    self.counter["Input files"]+=1
                    if fn.endswith(".gz"):
                        logger.debug("Opening file {} in gzip mode".format(fn))
                        fp = gzip.open(fn, "rt")
                    else:
                        logger.debug("Opening file {} in normal mode".format(fn))
                        fp = open(fn, "r")
                    f_list.append((fn,fp))

            return f_list

        else:
            raise SuperParserError ("Invalid file type")

class SuperParserError (Exception):
    """ Basic exception class for SuperParserError """
    pass
