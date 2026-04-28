#!/usr/bin/env python

import os
import re
import sys

from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import (
    _EMPTY,
    _EOL,
    _TAB
)

num = len


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s' % (message + _EOL * 3)

    stream.write(_EOL)
    stream.write('Usage: %s <orthogroups.tsv> [unplaced-regex]%s' % (
        os.path.basename(__file__), _EOL))
    stream.write(_EOL)
    stream.write('%s' % message)
    sys.exit(exitcode)
    

def main(argv):
    if len(argv) != 1 and len(argv) != 2:
        usage('Unexpected number of arguments')
    
    ortho = OrthoFinderOrthogroups(argv[0])
    if len(argv) == 2:
        unplaced = re.compile(argv[1])
    else:
        unplaced = re.compile('^$')
    
    M = dict()
    C = set()
    for g in range(num(ortho.groups)):
        for i in range(num(ortho.species)):
            for j in range(i, num(ortho.species)):
                for chr_i in ortho.groups[g][i]:
                    if chr_i is None:
                        continue
                    if unplaced.search(chr_i[len(ortho.species[i].id):]):
                        continue
                    for chr_j in ortho.groups[g][j]:
                        if chr_j is None:
                            continue
                        if unplaced.search(chr_j[len(ortho.species[j].id):]):
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

    C = sorted(C)  # dual purpose, sort and convert a list in one motion
    for i in range(num(ortho.species) - 1):
        C_i = list(filter(lambda c: c.startswith(ortho.species[i].id), C))
        for j in range(i + 1, num(ortho.species)):
            C_j = list(filter(lambda c: c.startswith(ortho.species[j].id), C))
            with open("%s-%s.matrix" % (ortho.species[i].id, ortho.species[j].id), 'wt') as ofile:
                ofile.write(_TAB.join(['Species'] + C_j) + _EOL)
                for chr_i in C_i:
                    ofile.write(chr_i)
                    for chr_j in C_j:
                        try:
                            count = M[chr_i][chr_j]
                        except:
                            count = 0
                            
                        ofile.write("%s%d" % (_TAB, count))

                    ofile.write(_EOL)

                
main(sys.argv[1:])
