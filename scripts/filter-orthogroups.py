#!/usr/bin/env python

import os
import sys
import getopt

from math import inf as _POS_INF
from og.core.io import is_stream
from og.core.common import map_loci_to_sequences
from og.core.parsers.config import SpeciesConfig
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr, is_placed
from og.constants import (
    _COMMENT,
    _EMPTY,
    _TAB,
    _EOL
)

_DEBUG = False

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'

_NEG_INF = -1.0 * _POS_INF
_STRICT = 1
_LENIENT = 2


num = len


    
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
    stream.write("  -N,--output-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then perform filtering.\n")
    stream.write("     Write sequence names to output (use `-n` for locus names).\n")
    stream.write("\n")
    stream.write("  -n,--map-to-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then perform filtering.\n")
    stream.write("     Write locus names to output (use `-N` for sequence names).\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized sequences to their\n")
    stream.write("     designated sequence names, not to their placed chromosome names.\n")
    stream.write("\n")
    stream.write("  -v,--invert-output\n")
    stream.write("     Output orthogroups failing the specified filters [passing]\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  1. The in.tsv file is an Orthogroups.tsv file with header defined.\n")
    stream.write("\n")
    # stream.write("  2. The newick-str argument is a Newick-formatted tree string that can be\n")
    # stream.write("     written to filter orthogroups by conditioning on the number of members\n")
    # stream.write("     (locus or sequence names) and species at each leaf node and internal\n")
    # stream.write("     node, respectively. In place of Newick branch lengths, however, the\n")
    # stream.write("     admissible number of members and species are specified using unsigned\n")
    # stream.write("     integer number ranges, consisting (inclusively) of the minimum number,\n")
    # stream.write("     a dash (`-`), then maximum number without any intervening whitespace.\n")
    # stream.write("     The minimum or maximum may be omitted to specify open ranges.\n")
    # stream.write("     For example, the tree below can be interpreted as follows:\n")
    # stream.write("         '((A:1, B:-4):1-2, (C:0-1, D):1-):2-4'\n")
    # stream.write("     Leaf node A must have one, and exactly one, member present; leaf node\n")
    # stream.write("     B may have up to four members (inclusive); the A+B subclade (here,\n")
    # stream.write("     represented as an internal node) requires one-to-two (inclusive)\n")
    # stream.write("     species to be be present. Leaf node C must be present at most once.\n")
    # stream.write("     The number of members on leaf node D is unrestricted; and subclade C+D\n")
    # stream.write("     requires at least one species be present. Because of the two subclade-\n")
    # stream.write("     specific constraints, the two-to-four species required at the root\n")
    # stream.write("     will always also be satisfied.\n")
    # stream.write("\n")    
    # #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    # #            0        10        20        30        40        50        60        70        80
    # stream.write("  3. A locus BED table is a two-column file specifying the species names\n")
    # stream.write("     (same as used in the Orthogroups.tsv header) and paths to locus BED\n")
    # stream.write("     files. The BED files must contain locus names in fourth column and the\n")
    # stream.write("     species names prepended to the sequence names (e.g., Hsa1 for chromosomes\n")
    # stream.write("     and HsaSca123 or HsaUn123 for unplaced scaffolds). If a locus BED table\n")
    # stream.write("     is given, then the number of sequences per species is counted as the\n")
    # stream.write("     members rather than locus names.\n")
    # stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    

def main(argv):
    short_flags = 'hiIuo:Nnv'
    long_flags = (
        'help',
        'output-file=',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'ignore-unlocalized',
        'output-sequence-names',
        'map-to-sequence-names',
        'invert-output'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    invert = False
    is_localized = is_placed
    map_seq_names = None
    output_seq_names = False
    output_file = sys.stdout
    ignore_unplaced = 0
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-u','--ignore-unlocalized'):
            is_localized = is_chr
        elif flag in ('-N','--output-sequence-names'):
            output_seq_names = map_seq_names = True
        elif flag in ('-n','--map-to-sequence-names'):
            map_seq_names = True
        elif flag in ('-v','--invert-output'):
            invert = True

    if not map_seq_names:
        output_seq_names = False
        ignore_unplaced = False

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    ortho = OrthoFinderOrthogroups(arguments[0])
    config = SpeciesConfig(arguments[1])
    tree = config.tree['ploidy']
    
    species_index = dict(zip(ortho.species, range(num(ortho.species))))

    for species in tree.terminal_nodes:
        if species.id not in species_index:
            raise KeyError("Species not found in orthologs file: '%s'" % (
                str(species.id)
            ))
    for species_id in ortho.species:
        if species_id not in config.species:
            raise KeyError("Species not found in YAML file: '%s'" % (
                str(species_id)
            ))

    output = OrthoFinderOrthogroups()
    output.species = ortho.species
    for group in range(num(ortho.groups)):
        counts = [0] * num(tree.nodes)
        if not map_seq_names:
            for i in tree.terminal_indices:
                species = tree.nodes[i]
                counts[i] = num(ortho.groups[group][species_index[species.id]])
        else:
            sequence_names = [tuple()] * num(ortho.species)
            for i in range(num(ortho.species)):
                sequence_names[i] = map_loci_to_sequences(
                    ortho.groups[group][i],
                    config.species[ortho.species[i]],
                    is_localized,
                    ignore_unplaced
                )
                if output_seq_names:
                    ortho.groups[group][i] = sorted(
                        sequence_names[i],
                        key=sequence_names[i].get,
                        reverse=True
                    )
                    
            for i in tree.terminal_indices:
                species = tree.nodes[i]
                sqnames = sequence_names[species_index[species.id]]
                counts[i] = sum(map(lambda c: int(c > 0), sqnames.values()))

                # The following conditions allow cells containing only unplaced
                # scaffolds to contribute toward presence count, and are
                # otherwise ignored in the presence of chromosome names.
                if ignore_unplaced == _LENIENT:
                    if counts[i] < 1 and \
                       sum(map(lambda c: int(c < 0), sqnames.values())) > 0:
                        counts[i] = species.length.maximum

        passes = True                        
        # anc = ancestor, dsc = descendant
        for anc_index, dsc_index in tree.get_edges(indices=True, reverse=True):
            dsc = tree.nodes[dsc_index]

            if dsc.length.minimum <= counts[dsc_index] <= dsc.length.maximum:
                if anc_index is not None:
                    counts[anc_index] += int(counts[dsc_index] > 0) \
                        if   (dsc.id in species_index) \
                        else counts[dsc_index]
            else:
                passes = False

        if invert:
            passes = not passes
                
        if passes:
            output.groups.append(ortho.groups[group])
            output.ids.append(ortho.ids[group])
        elif _DEBUG:
            sys.stderr.write(ortho.format_orthogroups_record(index=group) + _EOL)
            
    output.to_table(output_file)
    
    if is_stream(output_file):
        output_file.close()

        
main(sys.argv[1:])
