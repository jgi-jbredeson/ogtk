"""Module to provide IO support for reading and writing 
compressed and uncompressed files
"""

import io
import sys
import locale

__all__ = ('open','import_compression_module','STDIO')


from og.constants import _PYTHON_VERSION
from .constants import STDIO, COMPRESSION_NAME_MAP
from .filenames import infer_compression_format

# NOTES:
# - the compression modules (gzip, lzma, bz2file) accept file-like inputs,
#   even _io.TextIOWrapper objects for opening
#
# - the compression modules (gzip, lzma, bz2file) wrap the opened file
#   objects in _io.TextIOWrapper objects when 't' mode is requested
#
# - sys.stdin, sys.stderr, sys.stdout are _io.TextIOWrapper objects and
#   _io.TextIOWrapper cannot wrap another _io.TextIOWrapper object


def is_stream(fileobj):
    return hasattr(fileobj, '__next__') and hasattr(fileobj, 'close')


def infer_encoding(encoding=None):
    return locale.getpreferredencoding(False) if encoding is None else encoding



def import_compression_module(module):
    if type(module) == "module":
        # already loaded, pass-through
        module_object = module
    elif module in COMPRESSION_NAME_MAP:
        if COMPRESSION_NAME_MAP[module] == 'bgzip':
            import artisanal.core.io.bgzip as module_object
        elif COMPRESSION_NAME_MAP[module] == 'gzip':
            import gzip as module_object
        elif COMPRESSION_NAME_MAP[module] == 'bzip2':
            if _PYTHON_VERSION < (3,3):
                # stdlib bz2 module did not support multistream
                # prior to python 3.3; bz2file NOT in python stdlib
                import bz2file as module_object
            else:
                import bz2 as module_object
        elif COMPRESSION_NAME_MAP[module] == 'lzma':
            import lzma as module_object
        else:
            raise AssertionError("Should not be possible")
    elif module == 'io':
        import io as module_object
    else:
        raise ValueError("Unsupported compression format: %r" % module)
    return module_object



def open(filename, mode='rt', compresslevel=0, encoding=None, errors=None, newline=None, opener=None):
    """Open a gzip-/bgzip-/bzip2-/lzma-/xz-compressed file or uncompressed file
    in binary or text mode.

    The filename argument can be an actual filename (a str or bytes object), or
    an existing file object to read from or write to. Use "-" to open a file
    object to the appropriate stdin/stdout stream requested via mode.

    The mode argument can be "r", "rb", "w", "wb", "x", "xb", "a" or "ab" for
    binary mode, or "rt", "wt", "xt" or "at" for text mode. The default mode is
    "rt". The default compresslevel for compression streams is 9.

    For binary mode, the encoding, errors, and newline arguments must not be
    provided.

    For text mode, a compressed stream object is created, then wrapped in an
    io.TextIOWrapper instance with the specified encoding, error handling
    behavior, and line ending(s).

    To read from/write to compressed streams from stdin/stdout or existing
    file-like objects, the appropriate compression open() function object must
    be passed via the opener attribute (e.g., opener=gzip.open). When writing
    compresslevel must be set to a value between 1 and 9, inclusive (default
    is 9).

    This module relies on the io, gzip, lzma, bz2, bz2file, and pysam modules 
    internally, but other can be supplied via the opener attribute, which must
    return a callable object. bz2file is only required to for Python versions 
    earlier than 3.3.

    See Also: help(io.open)
    """
    compression = infer_compression_format(filename)
    encoding = infer_encoding(encoding)


    # if compression:
    #     if compresslevel < 1:
    #         compresslevel = 9
    #     modulename = COMPRESSION_NAME_MAP[compression]
    # if compresslevel:
    #     if not (1 <= compresslevel <= 9):
    #         raise ValueError("compresslevel must be between 1 and 9")        
    #     if compression is None and opener is None:
    #         # No file compression suffix, no opener; assume a stream
    #         raise ValueError("compresslevel specified for stream, but no opener given")
    # else:
    #     modulename = 'io'
    
    # Unless we explicitly request binary mode, default to
    # text mode. io, gzip, and lzma default to 'rt' mode, but
    # bz2file defaults to 'rb', so must set 'rt' explicitly
    # for a consistent interface
    binmode = 'b' in mode
    if 'r' in mode:
        mode = 'r'
    else:
        if 'a' in mode:
            mode = 'a'
        elif 'x' in mode:
            mode = 'x'
        elif 'w' in mode:
            mode = 'w'
        else:
            raise ValueError("Invalid mode: %r" % mode)

        if compression:
            if not compresslevel:
                compresslevel = 9
        if compresslevel:
            if not (1 <= compresslevel <= 9):
                raise ValueError("compresslevel must be between 1 and 9")
            if compression is None and opener is None:
                raise ValueError("compresslevel defined but no opener given")
        
    if binmode:
        mode += 'b'
    else:
        mode += 't'
        

    # Now define appropriate options for the opener:
    options = {
        'mode'   : mode,
        'errors' : errors,
        'newline': newline
    }
    if compresslevel:
        options['compresslevel'] = compresslevel
    elif not binmode:
        options['encoding'] = encoding

        
    # Lastly, open the file (if necessary)
    if filename == STDIO:
        stream = sys.stdin if 'r' in mode else sys.stdout
        if opener:
            # To use an opener() with our input/output stream,
            # use stream.buffer, as std streams decode bytes
            # to str under the hood automatically:
            # https://stackoverflow.com/questions/53245314/python-read-gzip-from-stdin
            return opener(stream.buffer, **options)
        else:
            # std uncompressed stream:
            return stream
        
    elif is_stream(filename):
        if opener is None:
            raise ValueError("file-like object defined but no opener given")
    else:
        opener = import_compression_module(
            COMPRESSION_NAME_MAP[compression] if compression else 'io'
        ).open
        
    return opener(filename, **options)



# def _load_module_object(module):
#     if type(module) == "module":
#         return module
#     else:
#         imports = {'module_object': None}
#         command = "import %s as module_object" % module
#         exec(command, imports)
#         return imports['module_object']
