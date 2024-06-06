#!/usr/bin/env python
# TODO: "Ignoring scaffolds" should allow single-scaffold membership

import os
import sys
import re
import getopt

from math import inf as _POS_INF
from og.core.io import is_stream
from og.core.parsers.newick import IntervalNewickTree
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import _COMMENT, _EMPTY, _TAB

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
    
def read_locus_bed(filename):
    records = dict()
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)
            
            assert num(fields) >= 4, \
                "Too few BED fields, expected at least four: '%s'" % line

            fields[0] = fields[0].strip()
            fields[3] = fields[3].strip()
            
            if fields[3] in records:
                sys.stderr.write("Duplicate locus ID: %s" % fields[3])
            else:
                records[fields[3]] = fields[0]
                
    return records



def read_locus_bed_table(filename):
    locus_bed = dict()
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(maxsplit=1)
            fields[0] = fields[0].strip()
            fields[1] = fields[1].strip()
            locus_bed[fields[0]] = read_locus_bed(fields[1])
            
    return locus_bed



def _map_loci_to_sequences(locus_list, locus_bed, species_id,
                           unplaced_re=None, ignore_unplaced=False):
    count = dict()
    for locus_name in locus_list:
        if locus_name not in locus_bed:
            raise KeyError("Locus ID not found in BED file: %s" % locus_name)

        sequence_name = locus_bed[locus_name]
        if sequence_name not in count:
            count[sequence_name] = 0        
        if unplaced_re and ignore_unplaced and \
           unplaced_re.search(sequence_name[len(species_id):]):
            count[sequence_name] -= 1
        else:
            count[sequence_name] += 1
    return count

    
    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <in.tsv> <newick-str>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -b,--locus-bed-table <file>\n")
    stream.write("     Table of species ID and path to a BED file for each species. Used to\n")
    stream.write("     map locus IDs in the input orthogroups file to sequence names.\n")
    stream.write("\n")
    stream.write("  -e,--regex-unplaced <regex>\n")
    stream.write("     Identify unplaced sequence using the specified regex. Takes effect\n")
    stream.write("     only when the `-b` option is also enabled [none]\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     Strictly ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains no chromosomal sequences, that cell\n")
    stream.write("     then contains no members. Takes effect only when the `-b` option is\n")
    stream.write("     also enabled.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     Leniently ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains only unplaced (ie, non-chomosomal)\n")
    stream.write("     sequences, that cell contains members. Takes effect only when the\n")
    stream.write("     `-b` option is also enabled.\n")
    stream.write("\n")
    stream.write("  -n,--map-to-sequence-names\n")
    stream.write("     Write sequence names to the output orthogroups table [locus IDs]\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
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
    stream.write("  2. The newick-str argument is a Newick-formatted tree string that can be\n")
    stream.write("     written to filter orthogroups by conditioning on the number of members\n")
    stream.write("     (locus or sequence IDs) and species at each leaf node and internal\n")
    stream.write("     node, respectively. In place of Newick branch lengths, however, the\n")
    stream.write("     admissible number of members and species are specified using unsigned\n")
    stream.write("     integer number ranges, consisting (inclusively) of the minimum number,\n")
    stream.write("     a dash (`-`), then maximum number without any intervening whitespace.\n")
    stream.write("     The minimum or maximum may be omitted to specify open ranges.\n")
    stream.write("     For example, the tree below can be interpreted as follows:\n")
    stream.write("         '((A:1, B:-4):1-2, (C:0-1, D):1-):2-4'\n")
    stream.write("     Leaf node A must have one, and exactly one, member present; leaf node\n")
    stream.write("     B may have up to four members (inclusive); the A+B subclade (here,\n")
    stream.write("     represented as an internal node) requires one-to-two (inclusive)\n")
    stream.write("     species to be be present. Leaf node C must be present at most once.\n")
    stream.write("     The number of members on leaf node D is unrestricted; and subclade C+D\n")
    stream.write("     requires at least one species be present. Because of the two subclade-\n")
    stream.write("     specific constraints, the two-to-four species required at the root\n")
    stream.write("     will always also be satisfied.\n")
    stream.write("\n")    
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  3. A locus BED table is a two-column file specifying the species IDs\n")
    stream.write("     (same as used in the Orthogroups.tsv header) and paths to locus BED\n")
    stream.write("     files. The BED files must contain locus IDs in fourth column and the\n")
    stream.write("     species IDs prepended to the sequence names (e.g., Hsa1 for chromosomes\n")
    stream.write("     and HsaSca123 or HsaUn123 for unplaced scaffolds). If a locus BED table\n")
    stream.write("     is given, then the number of sequences per species is counted as the\n")
    stream.write("     members rather than locus IDs.\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    

def main(argv):
    short_flags = 'hb:e:iIo:nv'
    long_flags = (
        'help',
        'output-file=',
        'regex-unplaced=',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'map-to-sequence-names',
        'locus-bed-table=',
        'invert-output'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    invert = False
    locus_table = None
    output_file = sys.stdout
    write_seq_names = False
    ignore_unplaced = 0
    regexp_unplaced = None  # re.compile('^(?:Sca|Un)', flags=re.IGNORECASE)
    for flag, value in options:
        if   flag in ('-h','--help'): usage(exitcode=0)
        elif flag in ('-o','--output-file'): output_file = value
        elif flag in ('-b','--locus-bed-table'): locus_table = value
        elif flag in ('-e','--regex-unplaced'): regexp_unplaced = re.compile(value)
        elif flag in ('-i','--ignore-unplaced-leniently'): ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'): ignore_unplaced = _STRICT
        elif flag in ('-n','--map-to-sequence-names'): write_seq_names = True
        elif flag in ('-v','--invert-output'): invert = True
        
    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    tree = IntervalNewickTree(arguments[1], default_length=(_NEG_INF, _POS_INF))
    ortho = OrthoFinderOrthogroups(arguments[0])

    species_index = dict(zip(ortho.species, range(num(ortho.species))))

    for species in tree.terminal_nodes:
        if species.id not in species_index:
            raise KeyError("Species not found in orthologs file: '%s'" % (
                str(species.id)
            ))

    if locus_table is None:
        write_seq_names = False
        ignore_unplaced = False
    else:
        locus_table = read_locus_bed_table(locus_table)
        for species_id in ortho.species:
            if species_id not in locus_table:
                raise KeyError("Species not found in locus BED file: '%s'" % (
                    str(species_id)
                ))
                    

    output = OrthoFinderOrthogroups()
    output.species = ortho.species
    for group in range(num(ortho.groups)):
        counts = [0] * num(tree.nodes)
        if locus_table is None:
            for i in tree.terminal_indices:
                species = tree.nodes[i]
                counts[i] = num(ortho.groups[group][species_index[species.id]])
        else:
            sequence_names = [tuple()] * num(ortho.species)
            for i in range(num(ortho.species)):
                sequence_names[i] = _map_loci_to_sequences(
                    ortho.groups[group][i],
                    locus_table[ortho.species[i]],
                    ortho.species[i],
                    regexp_unplaced,
                    ignore_unplaced
                )
                if write_seq_names:
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
            sys.stderr.write(ortho.format_orthogroups_record(index=group))
            
    output.to_table(output_file)
    
    if is_stream(output_file):
        output_file.close()

        
main(sys.argv[1:])
