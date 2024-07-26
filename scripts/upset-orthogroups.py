#!/usr/bin/env python

import os
import sys
import getopt

from og.core.io import open, STDIO
from og.core.common import _LENIENT, _STRICT
from og.core.common import map_loci_to_sequences
from og.core.common import filter_unplaced_sequences
from og.core.parsers.config import SpeciesConfig
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr, is_placed
from og.constants import _TAB, _SPACE, _EMPTY, _EOL

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'


num = len
_MEMBERSHIP = 0x1


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <in.tsv> <in.yaml>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -g,--group-by <enum>\n")
    stream.write("     Group output by 1 (membership) [0]\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     Strictly ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains no chromosomal sequences, that cell\n")
    stream.write("     then contains no members.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     Leniently ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains only unplaced (ie, non-chomosomal)\n")
    stream.write("     sequences, that cell contains members.\n")
    stream.write("\n")
    stream.write("  -M,--max-count <uint>\n")
    stream.write("     Count up to `-M` number of sequences per species and orthogroup.\n")
    stream.write("     Counts greater than this threshold are converted to the 'M' character\n")
    stream.write("     to denote 'Many' or 'Multiple' [2]\n")
    stream.write("\n")
    stream.write("  -n,--map-to-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then perform filtering.\n")
    stream.write("     Write locus names to output (use `-N` for sequence names).\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -p,--input-sequence-names\n")
    stream.write("     Locus names have already been mapped to their corresponding sequence\n")
    stream.write("     names in the input orthogroups file. Perform filtering accordingly.\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized sequences to their\n")
    stream.write("     designated sequence names, not to their placed chromosome names.\n")
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
    short_flags = 'hg:M:o:pniIu'
    long_flags = (
        'help',
        'group-by=',
        'max-count=',
        'output-file=',
        'input-sequence-names',
        'map-to-sequence-names',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'ignore-unlocalized'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    group_by = 0
    max_count = 2
    is_localized = is_placed
    map_seq_names = False
    input_seq_names = False
    output_file = STDIO
    ignore_unplaced = 0
    for flag, value in options:
        if flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-g','--group-by'):
            group_by = int(value)
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-M','--max-count'):
            max_count = int(value)
        elif flag in ('-n','--map-to-sequence-names'):
            map_seq_names = True
        elif flag in ('-p','--input-sequence-names'):
            input_seq_names = True
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-u','--ignore-unlocalized'):
            is_localized = is_chr

    if not (input_seq_names or map_seq_names):
        output_seq_names = False
        ignore_unplaced = False            
            
    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    ortho  = OrthoFinderOrthogroups(arguments[0])
    config = SpeciesConfig(arguments[1])
    output = open(output_file, 'w')
    
    for species_id in ortho.species:
        if species_id not in config.species:
            raise KeyError("Species not found in YAML file: '%s'" % (
                str(species_id)
            ))
        
    locus_counts = dict()
    ortho_counts = dict()
    num_species = num(ortho.species)
    for group in range(num(ortho.groups)):
        pattern = [0] * num_species
        sequence_names = list(ortho.groups[group])
        if map_seq_names:
            for s in range(num_species):
                sequence_names[s] = map_loci_to_sequences(
                    ortho.groups[group][s],
                    config.species[ortho.species[s]],
                    is_localized,
                    ignore_unplaced
                )
                 
        if input_seq_names or map_seq_names:
            for s in range(num_species):
                pattern[s] = num(filter_unplaced_sequences(
                    sequence_names[s],
                    config.species[ortho.species[s]],
                    is_localized,
                    ignore_unplaced
                ))
                pattern[s] = 'M' if pattern[s] > max_count else str(pattern[s])
        else:
            for s in range(num_species):
                pattern[s] = num(ortho.groups[group][s])
                pattern[s] = 'M' if pattern[s] > max_count else str(pattern[s])

        try:
            ortho_counts[tuple(pattern)] += 1
        except KeyError:
            ortho_counts[tuple(pattern)] = 1

    output.write(_SPACE.join(ortho.species) + _TAB + 'Count' + _EOL)
    if group_by & _MEMBERSHIP:
        group_patterns = dict()
        group_counts = dict()
        for pattern in ortho_counts:
            _pattern = tuple(sorted(pattern))
            if _pattern not in group_patterns:
                group_patterns[_pattern] = dict()
                group_counts[_pattern] = 0
            group_patterns[_pattern][pattern] = ortho_counts[pattern]
            group_counts[_pattern] += ortho_counts[pattern]

        for _pattern in sorted(group_patterns, key=group_counts.get, reverse=True):
            output.write("##total=%d\n" % group_counts[_pattern])
            for pattern in sorted(group_patterns[_pattern], key=group_patterns[_pattern].get, reverse=True):
                output.write(
                    _SPACE.join(pattern) + _TAB + str(ortho_counts[pattern]) + _EOL
                )
    else:
        for pattern in sorted(ortho_counts, key=ortho_counts.get, reverse=True):
            output.write(
                _SPACE.join(pattern) + _TAB + str(ortho_counts[pattern]) + _EOL
            )
            
    output.close()

    return 0
    

if __name__ == '__main__':
    main(sys.argv[1:])

