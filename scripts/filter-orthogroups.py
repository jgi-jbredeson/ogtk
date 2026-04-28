#!/usr/bin/env python

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'


import sys
import copy
import getopt

from math import inf as _POS_INF
from og.core.members import _LENIENT, _STRICT
from og.core.members import map_loci_to_sequences
from og.core.members import filter_unplaced_sequences
from og.core.compression import is_stream
from og.core.parsers.config import SampleConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
from og.constants import (
    _COMMENT,
    _EMPTY,
    _TAB,
    _EOL
)

_DEBUG = False
_NEG_INF = -1.0 * _POS_INF

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
    stream.write("  -m,--min-members <uint>\n")
    stream.write("     Minimum number of members permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -M,--max-members <uint>\n")
    stream.write("     Maximum number of members permitted per orthogroup [inf]\n")
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
    stream.write("  -p,--input-sequence-names\n")
    stream.write("     Locus names have already been mapped to their corresponding sequence\n")
    stream.write("     names in the input orthogroups file. Perform filtering accordingly.\n")
    stream.write("\n")
    stream.write("  -s,--min-samples <uint>\n")
    stream.write("     Minimum number of samples permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -S,--max-samples <uint>\n")
    stream.write("     Maximum number of samples permitted per orthogroup [inf]\n")
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
    # stream.write("     (locus or sequence names) and samples at each leaf node and internal\n")
    # stream.write("     node, respectively. In place of Newick branch lengths, however, the\n")
    # stream.write("     admissible number of members and samples are specified using unsigned\n")
    # stream.write("     integer number ranges, consisting (inclusively) of the minimum number,\n")
    # stream.write("     a dash (`-`), then maximum number without any intervening whitespace.\n")
    # stream.write("     The minimum or maximum may be omitted to specify open ranges.\n")
    # stream.write("     For example, the tree below can be interpreted as follows:\n")
    # stream.write("         '((A:1, B:-4):1-2, (C:0-1, D):1-):2-4'\n")
    # stream.write("     Leaf node A must have one, and exactly one, member present; leaf node\n")
    # stream.write("     B may have up to four members (inclusive); the A+B subclade (here,\n")
    # stream.write("     represented as an internal node) requires one-to-two (inclusive)\n")
    # stream.write("     samples to be be present. Leaf node C must be present at most once.\n")
    # stream.write("     The number of members on leaf node D is unrestricted; and subclade C+D\n")
    # stream.write("     requires at least one samples be present. Because of the two subclade-\n")
    # stream.write("     specific constraints, the two-to-four samples required at the root\n")
    # stream.write("     will always also be satisfied.\n")
    # stream.write("\n")    
    # #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    # #            0        10        20        30        40        50        60        70        80
    # stream.write("  3. A locus BED table is a two-column file specifying the samples names\n")
    # stream.write("     (same as used in the Orthogroups.tsv header) and paths to locus BED\n")
    # stream.write("     files. The BED files must contain locus names in fourth column and the\n")
    # stream.write("     samples names prepended to the sequence names (e.g., Hsa1 for chromosomes\n")
    # stream.write("     and HsaSca123 or HsaUn123 for unplaced scaffolds). If a locus BED table\n")
    # stream.write("     is given, then the number of sequences per samples is counted as the\n")
    # stream.write("     members rather than locus names.\n")
    # stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    

def main(argv):
    short_flags = 'ho:iIupm:M:Nns:S:v'
    long_flags = (
        'help',
        'output-file=',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'ignore-unlocalized',
        'input-sequence-names',
        'output-sequence-names',
        'map-to-sequence-names',
        'min-members=',
        'max-members=',
        'min-samples=','min-species=',
        'max-samples=','max-species=',
        'invert-output'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    invert = False
    is_placed = _placed
    min_members = 1
    max_members = _POS_INF
    min_samples = 1
    max_samples = _POS_INF
    map_seq_names = None
    input_seq_names = False
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
            is_placed = _localized
        elif flag in ('-m','--min-members'):
            min_members = int(value)
        elif flag in ('-M','--max-members'):
            max_members = int(value)
        elif flag in ('-N','--output-sequence-names'):
            output_seq_names = map_seq_names = True
        elif flag in ('-n','--map-to-sequence-names'):
            map_seq_names = True
        elif flag in ('-s','--min-samples','--min-species'):
            min_samples = int(value)
        elif flag in ('-S','--max-samples','--max-species'):
            max_samples = int(value)
        elif flag in ('-p','--input-sequence-names'):
            input_seq_names = True
        elif flag in ('-v','--invert-output'):
            invert = True

    if not (input_seq_names or map_seq_names):
        output_seq_names = False
        ignore_unplaced = False

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    ortho = OrthoFinderOrthogroups(arguments[0])
    config = SampleConfigFile(arguments[1], load_files=True, map_assigned_molecule=True)
    tree = config.tree['ploidy']
    
    sample_indices = {s.id:s.index for s in ortho.samples}

    for sample in tree.terminal_nodes:
        if sample.id not in sample_indices:
            raise KeyError(
                "Sample not found in orthologs file: '%s'" % str(sample.id)
            )
    for sample in ortho.samples:
        if sample.id not in config.samples:
            raise KeyError(
                "Sample not found in YAML file: '%s'" % str(sample.id)
            )

    output = OrthoFinderOrthogroups(samples=ortho.samples)
    for group in ortho.groups:
        counts = [0] * num(tree.nodes)
        sequence_names = list(group)
        if map_seq_names:
            for sample in ortho.samples:
                sequence_names[sample.index] = \
                    map_loci_to_sequences(
                        group[sample.index],
                        config.samples[sample.id],
                        is_placed,
                        ignore_unplaced
                    )
                if output_seq_names:
                    group[sample.index] = \
                        sequence_names[sample.index]

    
        for sample_index in tree.terminal_indices:
            sample = tree.nodes[sample_index]
            if input_seq_names or map_seq_names:
                counts[sample_index] = \
                    num(filter_unplaced_sequences(
                        sequence_names[sample_indices[sample.id]],
                        config.samples[sample.id],
                        is_placed,
                        ignore_unplaced
                    ))
            else:
                counts[sample_index] = \
                    num(group[sample_indices[sample.id]])

        passes = True
        if not (min_members <= sum(counts) <= max_members):
            passes = False
        if not (min_samples <= sum(map(bool, counts)) <= max_members):
            passes = False
            
        # anc = ancestor, dsc = descendant
        for anc_index, dsc_index in tree.get_edges(indices=True, reverse=True):
            dsc = tree.nodes[dsc_index]

            if dsc.length.minimum <= counts[dsc_index] <= dsc.length.maximum:
                if anc_index is not None:
                    counts[anc_index] += int(counts[dsc_index] > 0) \
                        if   (dsc.id in sample_indices) \
                        else counts[dsc_index]
            else:
                passes = False

        if invert:
            passes = not passes
                
        if passes:
            output.groups.append(group)
            
    output.to_file(output_file)
    
    if is_stream(output_file):
        output_file.close()

        
if __name__ == '__main__':
    main(sys.argv[1:])
