#!/usr/bin/env python

import os
import re
import sys

from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import (
    _EOL,
    _TAB
)

num = len
_UNPLACED = re.compile('(?:Sca|UN)')

def main(argv):
    ortho = OrthoFinderOrthogroups(argv[0])
    ofile = sys.stdout
    
    M = dict()
    C = set()
    for g in range(num(ortho.groups)):
        for i in range(num(ortho.species)):
            for j in range(i, num(ortho.species)):
                for chr_i in ortho.groups[g][i]:
                    if chr_i is None:
                        continue
                    if _UNPLACED.search(chr_i[len(ortho.species[i].id):]):
                        continue
                    for chr_j in ortho.groups[g][j]:
                        if chr_j is None:
                            continue
                        if _UNPLACED.search(chr_j[len(ortho.species[j].id):]):
                            continue
                        
                        if chr_i not in M:
                            M[chr_i] = {}
                        if chr_j not in M:
                            M[chr_j] = {}

                        if chr_j in M[chr_i]:
                            M[chr_i][chr_j] += 1
                        else:
                            M[chr_i][chr_j] = 1

                        if chr_i in M[chr_j]:
                            M[chr_j][chr_i] += 1
                        else:
                            M[chr_j][chr_i] = 1

                        C.add(chr_i)
                        C.add(chr_j)

    C = sorted(C)
    ofile.write(_TAB.join(['Species'] + C) + _EOL)
    for chr_i in C:
        ofile.write(chr_i)
        for chr_j in C:
            try:
                count = M[chr_i][chr_j]
            except:
                count = 0

            ofile.write("%s%d" % (_TAB, count))

        ofile.write(_EOL)
                
                
main(sys.argv[1:])
