
import os
import math

from og.core.generic import Enumerative
from og.core.tsv import read_to_list
from og.constants import _COMMA



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


