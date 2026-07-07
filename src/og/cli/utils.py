
import os
import math

from og.core.generic import Enumerative
from og.core.tsv import read_to_list
from og.constants import _COMMA,_SPACE

MAIN_ENTRY = False

def basename(filepath, suffix=None):
    filename = os.path.basename(filepath)
    if suffix and filename.endswith(suffix):
        filename = filename[:-len(suffix)]
    return filename


def get_exe(command, subcommand=None):
    if command:
        command = basename(command).replace('_','-')
    if subcommand:
        subcommand = basename(subcommand).replace('_','-')
    if MAIN_ENTRY:
        return _SPACE.join((command, subcommand))
    return command
    
    
def maxbits(maxsize=None):
    if maxsize is None:
        from sys import maxsize
    if maxsize < 0:
        maxsize = abs(maxsize)
    return int(math.log(maxsize, 2)+1)


def read_list(arg):
    items = []
    for item in arg.strip(_COMMA).split(_COMMA):
        if os.access(item, os.F_OK):  # a file
            items.extend(read_to_list(item))
        else:
            items.append(item)
    return items


