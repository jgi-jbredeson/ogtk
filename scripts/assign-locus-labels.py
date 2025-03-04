#!/usr/bin/env python

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Join labels to loci from OrthoFinder Orthogroups.tsv file'


import sys
import getopt

from math import inf as _POS_INF
from og.core.io import is_stream
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import (
    _COMMENT,
    _EMPTY,
    _TAB,
    _EOL
)

_DEBUG = False
_NEG_INF = -1.0 * _POS_INF
_STRAND = {'-': -1, '.': 0, '+': +1}

_STRAND_STR = '.+-'

num = len


def read_labels_file(labels_filename):
    labels = {}
    with open(labels_filename,'r') as labels_file:
        for line in labels_file:
            line = line.strip()
            if line == _EMPTY or \
               line.startswith('#'):
                continue

            fields = line.split()

            if num(fields) < 2:
                raise IndexError("Too few fields, at least two required")

            labels[fields[0]] = fields[1]

    return labels

    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <in.tsv> <labels.tsv> <species.id>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  1. The in.tsv file is an Orthogroups.tsv file with header defined.\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    

def main(argv):
    short_flags = 'ho:'
    long_flags = (
        'help',
        'output-file=',
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    output_file = sys.stdout
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-o','--output-file'):
            output_file = open(value, 'w')

    if num(arguments) != 3:
        usage('Unexpected number of arguments')

    ortho = OrthoFinderOrthogroups(arguments[0])
    labels = read_labels_file(arguments[1])
    species_id = arguments[2]
    
    
    if species_id not in ortho.species:
        raise KeyError("Species not found in Orthogroups file: '%s'" % (
            str(species_id)
        ))

    locus_labels = dict()
    species_index = ortho.species.index(species_id)
    for group in range(num(ortho.groups)):

        if ortho.ids[group] in labels:
            label = labels[ortho.ids[group]] 
            for locus_id in ortho.groups[group][species_index]:
                output_file.write("%s\t%s\n" % (locus_id, label))        
    
    if is_stream(output_file):
        output_file.close()

        
if __name__ == '__main__':
    main(sys.argv[1:])
