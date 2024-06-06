
from .constants import COMPRESSION_SUFFIX


def infer_compression_format(filename):
    """
    Infer from the filename extension the compression format used,
    or None.

    If passed a fileobj, returns None.
    """
    if isinstance(filename, str):
        tempname = filename.lower()
        for fileformat in sorted(COMPRESSION_SUFFIX):
            if tempname.endswith(COMPRESSION_SUFFIX[fileformat]):
                return fileformat
    return None



def strip_suffix(filename, suffixes=()):
    tempname = filename.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if tempname.endswith(suffix):
            filename = filename[:-len(suffix)]
    return filename


    
def split_suffix(filename, suffixes=()):
    tempname = filename.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if tempname.endswith(suffix):
            return filename[:-len(suffix)], filename[-len(suffix):]
    return filename, None



def strip_compression_suffix(filename):
    fileformat = infer_compression_format(filename)
    if fileformat is None:
        return filename
    else:
        return strip_suffix(filename, COMPRESSION_SUFFIX[fileformat])

