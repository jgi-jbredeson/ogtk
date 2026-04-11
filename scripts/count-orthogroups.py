#!/usr/bin/env python

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'


import sys
import getopt

from og.core.utils import count_items
from og.core.members import _LENIENT, _STRICT
from og.core.members import map_loci_to_sequences
from og.core.members import filter_unplaced_sequences
from og.core.compression import open, STDIO
from og.core.parsers.config import SpeciesConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
from og.constants import _TAB, _SPACE, _EMPTY, _EOL


num = len
_MEMBERSHIP = 0x1
_MULTIPLES = 0x2

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
    # stream.write("  -g,--group-by <enum>\n")
    # stream.write("     Group output by membership (=1), number of Multiples (=2) [0]\n")
    # stream.write("     (See the `--max-count` option below)\n")
    # stream.write("\n")
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
    stream.write("     Counts greater than this threshold are converted to the `M` character\n")
    stream.write("     to denote `Many` or `Multiple` [inf]\n")
    stream.write("\n")
    stream.write("  -n,--map-to-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then calculate upset.\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -p,--input-sequence-names\n")
    stream.write("     Locus names have already been mapped to their corresponding sequence\n")
    stream.write("     names in the input orthogroups file. Perform filtering accordingly.\n")
    stream.write("\n")
    stream.write("  -S,--count-species\n")
    stream.write("     Instead of member counts, write 1/0 for species presence/absence\n")
    stream.write("\n")
    stream.write("  -T,--output-totals\n")
    stream.write("     Output counts table with a column of row totals appended.\n")
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
    short_flags = 'hM:o:pSniITu'
    long_flags = (
        'help',
        'count-species',
        'max-count=',
        'output-file=',
        'output-totals',
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

    # group_by = 0
    max_count = float('inf')
    is_placed = _placed
    map_seq_names = False
    input_seq_names = False
    output_file = STDIO
    output_counts = True
    output_totals = False
    ignore_unplaced = 0
    for flag, value in options:
        if flag in ('-h','--help'):
            usage(exitcode=0)
        # elif flag in ('-g','--group-by'):
        #     group_by |= int(value)
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-T','--output-totals'):
            output_totals = True
        elif flag in ('-M','--max-count'):
            max_count = int(value)
        elif flag in ('-n','--map-to-sequence-names'):
            map_seq_names = True
        elif flag in ('-S','--count-species'):
            output_counts = False
        elif flag in ('-p','--input-sequence-names'):
            input_seq_names = True
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-u','--ignore-unlocalized'):
            is_placed = _localized

    if not (input_seq_names or map_seq_names):
        output_seq_names = False
        ignore_unplaced = False            
            
    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    ortho  = OrthoFinderOrthogroups(arguments[0])
    config = SpeciesConfigFile(arguments[1], load_files=True, map_assigned_molecule=True)
    output = open(output_file, 'w')
    
    for species_id in ortho.species:
        if species_id not in config.species:
            raise KeyError("Species not found in YAML file: '%s'" % (
                str(species_id)
            ))

    totals_label = ['Total'] if output_totals else []
    output.write('\t'.join(['Orthogroup'] + ortho.species + totals_label) + '\n')
    member_counts = dict()
    ortho_counts = dict()
    num_species = num(ortho.species)
    for group_i in range(num(ortho.groups)):
        num_members = 0
        pattern = [0] * num_species
        sequence_names = [None] * num_species
        if map_seq_names:
            for species_i in range(num_species):
                sequence_names[species_i] = map_loci_to_sequences(
                    ortho.groups[group_i][species_i],
                    config.species[ortho.species[species_i]],
                    is_placed,
                    ignore_unplaced=False
                )
        else:  # either loci or pre-mapped sequences:
            for species_i in range(num_species):
                sequence_names[species_i] = count_items(
                    ortho.groups[group_i][species_i]
                )
                 
        if input_seq_names or map_seq_names:
            for species_i in range(num_species):
                names = filter_unplaced_sequences(
                    sequence_names[species_i],
                    config.species[ortho.species[species_i]],
                    is_placed,
                    ignore_unplaced,
                    aggregate_unplaced=False
                )
                if output_counts:
                    pattern[species_i] = 'M' if num(names) > max_count else num(names)
                    num_members += num(names)
                else:
                    pattern[species_i] = int(bool(names))
                    num_members += pattern[species_i]
        else:
            for species_i in range(num_species):
                names = ortho.groups[group_i][species_i]
                if output_counts:
                    pattern[species_i] = 'M' if num(names) > max_count else num(names)
                    num_members += num(names)
                else:
                    pattern[species_i] = int(bool(names))
                    num_members += pattern[species_i]
        
        output.write(
            '%s\t%s%s\n' % (
                ortho.ids[group_i],
                '\t'.join(map(str, pattern)),
                '\t%d' % (num_members) if output_totals else ''
            )
        )
            
    output.close()

    return 0
    

if __name__ == '__main__':
    main(sys.argv[1:])

