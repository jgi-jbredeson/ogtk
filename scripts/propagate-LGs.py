#!/usr/bin/env python

import os
import sys

from og.core.io import open

def fill_missing_cns(groups, letters=[]):
    if len(letters) < 1:
        return groups
    if len(letters) > 1:
        return groups
    letter = list(letters)[0]
    for group in groups:
        group[-1] = letter
    return groups


def main(argv):
    orthofile = open(argv[0],'r')
    o = 0  # column offset
    groups = list()
    letters = set()
    prev_og = None
    for line in orthofile:
        line = line.lstrip().rstrip('\r\n')

        if line == '' or \
           line.startswith('#'):
            continue

        if line.startswith('Cluster') or \
           line.startswith('Orthogroup'):
            print(line)
            o = 0
            continue
        elif line.startswith('Count'):
            print(line)
            o = 1
            continue

        fields = line.split('\t')

        curr_og = fields[o].split('.')[0]
        if prev_og is not None:
            if curr_og != prev_og:
                for group in fill_missing_cns(groups, letters):
                    print('\t'.join(group))
                letters = set()
                groups = list()

        fields[-1] = fields[-1].strip()
        if len(fields[-1]) > 0:
            letters.add(fields[-1])
        groups.append(fields)
        
        prev_og = curr_og

    if groups:
        for group in fill_missing_cns(groups, letters):
            print('\t'.join(group))
        letters = set()
        groups = list()

        
main(sys.argv[1:])
