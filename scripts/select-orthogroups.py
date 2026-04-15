#!/usr/bin/env python

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Select OrthoFinder Orthogroups.tsv file members'


import sys
import getopt
from og.core.utils import index_list
from og.core.parsers.tsv import read_to_list
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import _EMPTY


num = len



def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <in.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -L,--list-file <file>\n")
    stream.write("     Input file listing orthogroup members to include\n")
    stream.write("\n")
    stream.write("  -m,--min-members <int>\n")
    stream.write("     Minimum number of members per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -s,--min-species <int>\n")
    stream.write("     Minimum number of species per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -r,--no-remove-members\n")
    stream.write("     Do not remove members, only exhausted orthogroups\n")
    stream.write("\n")
    stream.write("  -v,--invert\n")
    stream.write("     Write lines that exhaust their members\n")
    stream.write("\n")
    stream.write("  -x,--exclude\n")
    stream.write("     Exclude members listed in list file\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    stream.write("\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n%s" % message)
    sys.exit(exitcode)
    
    
def main(argv):
    short_flags = 'hL:m:s:rxv'
    long_flags = (
        'help',
        'list-file=',
        'min-members=',
        'min-species=',
        'no-remove-members',
        'exclude',
        'invert'
    )

    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    invert = False
    exclude = False
    keepset = []
    min_members = 1
    min_species = 1
    remove_members = True
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-L','--list-file'):
            keepset = set(read_to_list(value))
        elif flag in ('-m','--min-members'):
            min_members = value
        elif flag in ('-s','--min-species'):
            min_species = value
        elif flag in ('-r','--no-remove-members'):
            remove_members = False
        elif flag in ('-v','--invert'):
            invert = True
        elif flag in ('-x','--exclude'):
            exclude = True
    try:
        min_members = int(min_members)
    except:
        usage("--min-members must be a positive integer")
    try:
        min_species = int(min_species)
    except:
        usage("--min-species must be a positive integer")
    if min_members < 1:
        usage("--min-members must be a positive integer")
    if min_species < 1:
        usage("--min-species must be a positive integer")
        
    if num(arguments) == 0:
        usage()
    if num(arguments) != 1:
        usage('Unexpected number of arguments')
    
    ortho = OrthoFinderOrthogroups(arguments[0])
    
    speciesA_index = index_list(ortho.species)
    groupsA_index = index_list(ortho.ids)

    output = OrthoFinderOrthogroups()
    output.species = ortho.species

    for g in range(num(ortho.groups)):
        members_counts = [0] * num(ortho.species)
        species_counts = [0] * num(ortho.species)

        group = []
        for s in range(num(ortho.species)):
            members = []
            if ortho.groups[g][s]:
                for member in ortho.groups[g][s]:
                    present = member in keepset
                    if exclude:
                        present = not present
                    if present:
                        members.append(member)
                    present = int(present)
                    members_counts[s] += present
                    species_counts[s] = present
            group.append(members)

        passes = ((sum(species_counts) >= min_species) and
                  (sum(members_counts) >= min_members))
        if invert:
            passes = not passes
        if passes:
            output.groups.append(group if remove_members else ortho.groups[g])
            output.ids.append(ortho.ids[g])
            
    output.to_table(sys.stdout)


if __name__ == '__main__':
    main(sys.argv[1:])

