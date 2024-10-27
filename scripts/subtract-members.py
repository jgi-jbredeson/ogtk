#!/usr/bin/env python

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'


import sys
from og.core.utils import index_list
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import _EMPTY


num = len


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else "\nERROR: %s\n\n\n" % message
    stream.write("Usage: %s <A.tsv> <B.tsv>\n%s" % (__program__, message))
    sys.exit(exitcode)

    
def main(argv):
    if num(argv) == 0:
        usage()
    if num(argv) != 2:
        usage('Unexpected number of arguments')
    
    orthoA = OrthoFinderOrthogroups(argv[0])
    orthoB = OrthoFinderOrthogroups(argv[1])

    speciesA_index = index_list(orthoA.species)
    speciesB_index = index_list(orthoB.species)

    groupsA_index = index_list(orthoA.ids)
    groupsB_index = index_list(orthoB.ids)

    for g in range(num(orthoA.groups)):
        idA = orthoA.ids[g]
        if idA not in groupsB_index:
            continue
        for s in range(num(orthoA.species)):
            speciesA = orthoA.species[s]
            if speciesA not in speciesB_index:
                continue
            orthoA.groups[g][s] = tuple(set(orthoA.groups[g][s]) \
                - set(orthoB.groups[groupsB_index[idA]][speciesB_index[speciesA]]))

    orthoA.to_table(sys.stdout)


if __name__ == '__main__':
    main(sys.argv[1:])

