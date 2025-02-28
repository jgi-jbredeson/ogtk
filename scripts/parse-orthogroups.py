#!/usr/bin/env python

import os
import sys

from getopt import getopt, GetoptError
#from fastx.utils.io import load_listfile    
from og.core.io import is_stream
from og.core.utils import index_list
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import (
    _COMMA,
    _EMPTY,
    _SPACE,
    _EOL,
    _TAB
)

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Manipulate OrthoFinder orthogroups files'

_LF = '\n'
_CR = '\r'

num = len


def lineparser(line):
    return line.rstrip(_EOL).split(_TAB)[0]


def load_listfile(listfiles, parser=lambda x: x):
    """Load a non-redundant list of items from listfiles.

    When multiple instances of the same value item are encountered,
    just the first is kept.

    parser is a callable that inputs a single string and returns 
    an object to be appended to the end of the growing list.
    """
    if isinstance(listfiles, str):
        listfiles = (listfiles,)
    itemset = set()
    itemlist = list()
    for listfile in listfiles:
        inputfile = open(listfile, 'r')
        for line in inputfile:
            item = parser(line)

            if item in itemset:
                continue
        
            itemlist.append(item)
            itemset.add(item)

        inputfile.close()

    return itemlist


def read_order_file(filename):
    return load_listfile(filename, parser=lineparser)
    # if "Orthogroup" not in order:
    #     order = ["Orthogroup"] + order
    # return order


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <Orthogroups.tsv>\n" % (
        os.path.basename(sys.argv[0])
    ))
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -d,--prefix-delim <char>\n")
    stream.write("     When prefixing species names to locus IDs, seperate them using char [|]\n")
    stream.write("\n")
    stream.write("  -D,--replace-delim <char>\n")
    stream.write("     To make IDs safe, replace existing specified char in locus IDs prior to\n")
    stream.write("     prefixing species ID [|]\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -O,--output-type <str>\n")
    stream.write("     Output orthogroups in specified format. [F], OrthoFinder; V, OrthoVenn\n")
    stream.write("\n")
    stream.write("  -p,--prefix-species-names\n")
    stream.write("     Prepend the locus IDs with the species names declared in the header.\n")
    stream.write("\n")
    stream.write("  -S,--species-order-file <file>\n")
    stream.write("     Input file listing (one per line) the desired output species order.\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this usage message\n")
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  - ClusterVenn3 requires `.txt` file suffix\n")
    stream.write("\n\n%s" % message)
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    sys.exit(exitcode)
    
    
def main(argv):
    err = sys.stderr
    output_file = sys.stdout
    species_delim = '|'
    replace_delim = species_delim
    species_prefix = False
    oldspecies = None
    oldindices = None
    newspecies = None
    newindices = None
    output_type = 'F'
    valid_output_types = set('FV')
    long_flags = (
        'outfile=',
        'output-file=',
        'output-type=',
        'species-order=',
        'species-order-file=',
        'replace-delim',
        'prefix-delim=',
        'prefix-species-names',
        'orthovenn',
        'help'
    )
    short_flags = 'd:o:O:D:pS:vh'
    
    try:
        options, arguments = getopt(argv, short_flags, long_flags)
    except GetoptError as message:
        usage(message)

    for flag, value in options:
        if   flag in ('-o','--outfile','--output-file'):
            output_file = open(value, 'w')
        elif flag in ('-O','--output-type'):
            output_type = value
        elif flag in ('-d','--prefix-delim'):
            species_delim = value
        elif flag in ('-D','--replace-delim'):
            replace_delim = value
        elif flag in ('-p','--prefix-species-names'):
            species_prefix = True
        elif flag in ('-S','--species-order','--species-order-file'):
            newspecies = read_order_file(value)
            newindices = index_list(newspecies)
        elif flag in ('-v','--orthovenn'):
            output_type = 'V'
        elif flag in ('-h','--help'):
            usage(exitcode=0)

    if output_type not in valid_output_types:
        usage("Unsupported output type: `%s`" % output_type)
            
    if len(arguments) != 1:
        usage("Unexpected number of arguments")

    #err.write("Outputting species in the following order:\n")
    #err.write("  %s\n" % ', '.join(newspecies))
        
    ortho = OrthoFinderOrthogroups(arguments[0])
    
    oldspecies = ortho.species
    oldindices = index_list(oldspecies)

    if newspecies is None:
        newspecies = oldspecies
        newindices = oldindices
    else:
        for species in newspecies:
            if species not in oldindices:
                raise KeyError("Species not found in orthogroups file: `%s`" % species)
            
    ortho.species = newspecies    
    for g in range(num(ortho.groups)):
        oldgroup = ortho.groups[g]
        newgroup = [None] * num(newspecies)
        for species in newspecies:
            newgroup[newindices[species]] = oldgroup[oldindices[species]]
            if newgroup[newindices[species]] is None:
                newgroup[newindices[species]] = tuple()
        ortho.groups[g] = newgroup

    if species_prefix:
        for g in range(num(ortho.groups)):
            for s in range(num(ortho.species)):
                ortho.groups[g][s] = tuple(
                    map(lambda L: ortho.species[s] + species_delim + L.replace(species_delim,replace_delim),
                        ortho.groups[g][s])
                )
        
    if output_type == 'V':
        orthogroups = []
        for g in range(num(ortho.groups)):
            n = 0
            for s in range(num(ortho.species)):
                n += int(len(ortho.groups[g][s]) > 0)
            orthogroups.append((-n, ortho.ids[g], ortho.groups[g]))

        orthogroups.sort()

        # output_file.write(_TAB.join(ortho.species) + _EOL)
        for g in range(num(orthogroups)):
            ortho.ids[g] = orthogroups[g][1]
            ortho.groups[g] = orthogroups[g][2]
            sep = _EMPTY
            for s in range(num(ortho.species)):
                if len(ortho.groups[g][s]) > 0:
                    output_file.write("%s%s" % (sep, _TAB.join(ortho.groups[g][s])))
                    sep = _TAB
            output_file.write(_LF)
    else:
        ortho.to_table(output_file)

    if is_stream(output_file):
        output_file.close()
        
main(sys.argv[1:])
