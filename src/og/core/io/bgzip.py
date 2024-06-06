"""Module for reading and writing bgzip-compressed files.

"""

# This code copied with modification from:
# https://github.com/nvawda/bz2file/blob/master/bz2file.py

__all__ = ('BGZFile','open')

import io
import sys

from .constants import _EMPTY
from pysam import BGZFile as _BGZFile

_STR_TYPES = (str, unicode) if (str is bytes) else (str, bytes)


def BGZFile(filename, mode="r", buffering=None, compresslevel=9):
    """Open a bgzip-compressed file.

    If filename is a str, bytes or unicode object, it gives the name
    of the file to be opened. Otherwise, it should be a file object,
    which will be used to read or write the compressed data.
        
    mode can be 'r' for reading (default), 'w' for (over)writing,
    'x' for creating exclusively, or 'a' for appending. These can
    equivalently be given as 'rb', 'wb', 'xb', and 'ab'.
    
    buffering is ignored. Its use is deprecated.
        
    compresslevel is currently ignored, but may be exposed to the user
    in the future. If mode is 'w', 'x' or 'a', compresslevel can be a 
    number between 1 and 9 specifying the level of compression: 
    1 produces the least compression, and 9 (default) produces the most
    compression.

    If mode is 'r', the input file may be the concatenation of
    multiple compressed streams.
    """

    if not (1 <= compresslevel <= 9):
        raise ValueError("compresslevel must be between 1 and 9")

    if mode in (_EMPTY, "r", "rb"):
        mode = "rb"
    elif mode in ("w", "wb"):
        mode = "wb"
    elif mode in ("x", "xb"):
        mode = "xb"
    elif mode in ("a", "ab"):
        mode = "ab"
    else:
        raise TypeError("filename must be a %s or %s object, or a file" %
                        (_STR_TYPES[0].__name__, _STR_TYPES[1].__name__))

    return _BGZFile(filename, mode)


def open(filename, mode="rb", compresslevel=9,
         encoding=None, errors=None, newline=None):
    """Open a bgzip-compressed file in binary or text mode.

    The filename argument can be an actual filename (a str, bytes or unicode
    object), or an existing file object to read from or write to.

    The mode argument can be "r", "rb", "w", "wb", "x", "xb", "a" or
    "ab" for binary mode, or "rt", "wt", "xt" or "at" for text mode.
    The default mode is "rb", and the default compresslevel is 9.

    For binary mode, this function is equivalent to the BGZFile
    constructor: BGZFile(filename, mode, compresslevel). In this case,
    the encoding, errors and newline arguments must not be provided.
    
    For text mode, a BGZFile object is created, and wrapped in an
    io.TextIOWrapper instance with the specified encoding, error
    handling behavior, and line ending(s).
    """
    if "t" in mode:
        if "b" in mode:
            raise ValueError("Invalid mode: %r" % (mode,))
    else:
        if encoding is not None:
            raise ValueError("Argument 'encoding' not supported in binary mode")
        if errors is not None:
            raise ValueError("Argument 'errors' not supported in binary mode")
        if newline is not None:
            raise ValueError("Argument 'newline' not supported in binary mode")

    bz_mode = mode.replace("t", _EMPTY)
    binary_file = BGZFile(filename, bz_mode, compresslevel=compresslevel)

    if "t" in mode:
        return io.TextIOWrapper(binary_file, encoding, errors, newline)
    else:
        return binary_file

    
