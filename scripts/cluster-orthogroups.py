#!/usr/bin/env python3
# using python3 ensures unlimited max integer (>> 2**63 - 1)

# TODO: write up phylogenetic grouping algorithm:
#       - single-chromosome mismatches (fissions/fusions)
#       - consistent (sub-)clade patterns (= translocations; will tend to have few counts)
#         e.g., fish vs bony vert specific chromosome mismatches
# TODO: merge unplaced+singletons group with reference groups, allowing scaffolds to contribute to count


import os
import sys
import re
import getopt

from math import inf as _POS_INF
from og.core.members import _LENIENT, _STRICT
from og.core.members import map_loci_to_sequence_counts
from og.core.members import filter_unplaced_sequences
from og.core.compression import is_stream
from og.core.parsers.config import SampleConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed

from og.constants import (
    _COMMA,
    _COMMENT,
    _EMPTY,
    _SPACE,
    _TAB
)

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'


_CLUSTER_ID = "CL{0:05d}".format
_PATTERN_ID = "{0:s}.{1:d}".format
_COMMASPACE = _COMMA + _SPACE
_EPSILON = 1e-6
_NEG_INF = -1.0 * _POS_INF

_STDOUT = sys.stdout
_STDERR = sys.stderr


num = len



def _as_tuples(listobj):
    return tuple(map(lambda l: tuple(l) if l else None, listobj))


def _as_set(listobj):
    _set = set()
    for item in listobj:
        if item is not None:
            _set.update(item)
    return _set

                                               
def calc_dist(list_u, list_v, countgaps=False):
    matches = 0.0
    mismatches = 0.0
    gaps = 0
    for i in range(len(list_u)):
        len_u = len(list_u[i])
        len_v = len(list_v[i])
        if len_u > 0 and len_v > 0:
            if list_u[i] == list_v[i]:
                matches += 1.0
            else:
                num_intersect = num(list_u[i].intersection(list_v[i]))
                if num_intersect > 0:
                    matches += num_intersect / num(list_u[i].union(list_v[i]))
                else:
                    mismatches += 1.0
        elif len_u == 0 and len_v == 0:
            continue
        elif countgaps:
            gaps += 1
    return float(mismatches + gaps) / float(max(1, matches + mismatches + gaps))
                

def _notNone(obj):
    return obj is not None




def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <in.tsv> <in.conf>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    # stream.write("  -b,--locus-bed-table <file>\n")
    # stream.write("     Table of sample ID and BED path\n")
    stream.write("\n")
    stream.write("  -c,--output-cluster-counts-file <file>\n")
    stream.write("     Write distinct cluster patterns with counts to file.\n")
    stream.write("\n")
    stream.write("  -C,--output-cluster-counts-file-all <file>\n")
    stream.write("     Write all distinct cluster patterns with counts to file.\n")
    stream.write("\n")
    # stream.write("  -d,--min-distance <ufloat>\n")
    # stream.write("     Minimum distance between orthogroups to cluster [0.0]\n")
    # stream.write("\n")
    # stream.write("  -D,--max-distance <ufloat>\n")
    # stream.write("     Maximum distance between orthogroups to cluster [1.0]\n")
    # stream.write("\n")
    # stream.write("  -e,--regex-unplaced <regex>\n")
    # stream.write("     Identify unplaced sequence using the specified regex. Takes effect\n")
    # stream.write("     only when the `-b` option is also enabled [none]\n")
    # stream.write("\n")
    stream.write("  -F,--output-cluster-map-file <file>\n")
    stream.write("     Write each cluster ID and member orthogroup ID to file.\n")
    stream.write("\n")    
    stream.write("  -g,--min-orthogroups <uint>\n")
    stream.write("     Minimum number of orthogroups per cluster to output [1]\n")
    stream.write("\n")
    stream.write("  -G,--max-orthogroups <uint>\n")
    stream.write("     Maximum number of orthogroups per cluster to output [inf]\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     Strictly ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains no chromosomal sequences, that cell\n")
    stream.write("     then contains no members.\n")  # Takes effect only when the `-b` option is\n")
    # stream.write("     also enabled.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     Leniently ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains only unplaced (ie, non-chomosomal)\n")
    stream.write("     sequences, that cell contains members.\n")  # Takes effect only when the\n")
    # stream.write("     `-b` option is also enabled.\n")
    stream.write("\n")
    stream.write("  -m,--min-members <uint>\n")
    stream.write("     Minimum number of overlapping members between orthogroup sets [2]\n")
    stream.write("\n")
    stream.write("  -M,--max-members <uint>\n")
    stream.write("     Maximum number of overlapping members between orthogroup sets [inf]\n")
    stream.write("\n")
    stream.write("  -N,--output-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then perform clustering.\n")
    stream.write("     Write sequence names to output (use `-n` for locus names).\n")
    stream.write("\n")
    stream.write("  -n,--map-to-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then perform clustering.\n")
    stream.write("     Write locus names to output (use `-N` for sequence names).\n")
    stream.write("\n")    
    stream.write("  -s,--min-samples <uint>\n")
    stream.write("     Minimum number of samples permitted per orthogroup [2]\n")
    stream.write("\n")
    stream.write("  -S,--max-samples <uint>\n")
    stream.write("     Maximum number of samples permitted per orthogroup [inf]\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized sequences to their\n")
    stream.write("     designated sequence names, not to their placed chromosome names.\n")
    stream.write("\n")    
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit.\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  1. The in.tsv file is an Orthogroups.tsv file with header defined and\n")
    stream.write("     assumes each orthogroup contains sequence names as members or that the\n")
    stream.write("     `-b` option is also enabled to map input locus names to sequence names.\n")
    stream.write("\n")
    # stream.write("  2. The newick-str argument is a Newick-formatted tree string that can be\n")
    # stream.write("     written to filter orthogroups by conditioning on the number of members\n")
    # stream.write("     (locus or sequence IDs) and sample at each leaf node and internal\n")
    # stream.write("     node, respectively. In place of Newick branch lengths, however, the\n")
    # stream.write("     admissible number of members and sample are specified using unsigned\n")
    # stream.write("     integer number ranges, consisting (inclusively) of the minimum number,\n")
    # stream.write("     a dash (`-`), then maximum number without any intervening whitespace.\n")
    # stream.write("     The minimum or maximum may be omitted to specify open ranges.\n")
    # stream.write("     For example, the tree below can be interpreted as follows:\n")
    # stream.write("         '((A:1, B:-4):1-2, (C:0-1, D):1-):2-4'\n")
    # stream.write("     Leaf node A must have one, and exactly one, member present; leaf node\n")
    # stream.write("     B may have up to four members (inclusive); the A+B subclade (here,\n")
    # stream.write("     represented as an internal node) requires one-to-two (inclusive)\n")
    # stream.write("     sample to be be present. Leaf node C must be present at most once.\n")
    # stream.write("     The number of members on leaf node D is unrestricted; and subclade C+D\n")
    # stream.write("     requires at least one sample be present. Because of the two subclade-\n")
    # stream.write("     specific constraints, the two-to-four sample required at the root\n")
    # stream.write("     will always also be satisfied.\n")
    # #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    # #            0        10        20        30        40        50        60        70        80
    # stream.write("\n")
    # stream.write("  3. A locus BED table is a two-column file specifying the sample IDs\n")
    # stream.write("     (same as used in the Orthogroups.tsv header) and paths to locus BED\n")
    # stream.write("     files. The BED files must contain locus IDs in fourth column and the\n")
    # stream.write("     sample IDs prepended to the sequence names (e.g., Hsa1 for chromosomes\n")
    # stream.write("     and HsaSca123 or HsaUn123 for unplaced scaffolds). If a locus BED table\n")
    # stream.write("     is given, then the number of sequences per sample is counted as the\n")
    # stream.write("     members rather than locus IDs.\n")
    # stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    
def main(argv):
    short_flags = 'hc:C:d:D:F:g:G:iIm:M:nNs:S:u'
    long_flags = (
        'help',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'ignore-unlocalized',
        'min-members=',
        'max-members=',
        'min-orthogroups',
        'max-orthogroups',
        'min-samples=','min-species=',
        'max-samples=','max-species=',
        'min-distance=',
        'max-distance=',
        'output-cluster-counts-file-all=',
        'output-cluster-counts-file=',
        'output-cluster-map-file=',
        'output-sequence-names',
        'map-to-sequence-names'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    min_dist = 0.0
    max_dist = 1.0
    min_members = 2
    max_members = _POS_INF
    min_samples = 2
    max_samples = _POS_INF
    min_orthogroups = 1
    max_orthogroups = _POS_INF
    cluster_map_file = False
    cluster_counts_all_file = False
    cluster_counts_mrg_file = False
    is_placed = _placed
    ignore_unplaced = False
    map_seq_names = False
    output_seq_names = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-u','--ignore-unlocalized'):
            is_placed = _localized
        elif flag in ('-m','--min-members'):
            min_members = int(float(value))
        elif flag in ('-M','--max-members'):
            max_members = int(float(value))
        elif flag in ('-g','--min-orthogroups'):
            min_orthogroups = int(float(value))
        elif flag in ('-G','--max-orthogroups'):
            max_orthogroups = int(float(value))
        elif flag in ('-s','--min-samples','--min-species'):
            min_samples = int(float(value))
        elif flag in ('-S','--max-samples','--max-species'):
            max_samples = int(float(value))
        elif flag in ('-d','--min-distance'):
            min_dist = float(value)
        elif flag in ('-D','--max-distance'):
            max_dist = float(value)
        elif flag in ('-n','--map-to-sequence-names'):
            map_seq_names = True
        elif flag in ('-N','--output-sequence-names'):
            output_seq_names = map_seq_names = True
        elif flag in ('-C','--output-cluster-counts-file-all'):
            cluster_counts_all_file = open(value, 'wt')
        elif flag in ('-c','--output-cluster-counts-file'):
            cluster_counts_mrg_file = open(value, 'wt');
        elif flag in ('-F','--output-cluster-map-file'):
            cluster_map_file = open(value, 'wt')

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    loc_ortho  = OrthoFinderOrthogroups(arguments[0])
    config = SampleConfigFile(arguments[1], load_files=True, map_assigned_molecule=True)
    clust  = config.tree['ploidy']
    ofile  = sys.stdout

    num_samples = num(loc_ortho.samples)
    samples_index = {s.id:s.index for s in loc_ortho.samples}

    if min_dist >= 1.0:
        min_dist = min_dist / num_samples - _EPSILON
    if max_dist >= 1.0:
        max_dist = max_dist / num_samples + _EPSILON

    if map_seq_names:
        for sample in loc_ortho.samples:
            if sample.id not in config.samples:
                raise KeyError(
                    "Sample not found in conf file: '%s'" % str(sample.id)
                )

        seq_ortho = OrthoFinderOrthogroups(samples=loc_ortho.samples)
        # ortho.samples = loc_ortho.samples
        # ortho.ids = loc_ortho.ids
        # ortho.groups = [None] * num(loc_ortho.groups)
        for loc_group in loc_ortho.groups:
            # ortho.groups[group] = [None] * num_samples
            seq_group = seq_ortho.new_group(append=True)
            for sample in loc_ortho.samples:
                seq_group[sample.index] = \
                    list(map_loci_to_sequence_counts(
                        loc_group[sample.index],
                        config.samples[sample.id],
                        is_placed,
                        ignore_unplaced
                    ))
        ortho = seq_ortho
    else:
        ortho = loc_ortho
    # TODO:
    #  Instead of using strings in sets, index gene-containing sequences
    #  into bit arrays and use std set operations. Must use numpy.array,
    #  as set operations are not supported by built-in list objects.
    #
    #  Store bit masks for all genome chr indices
    #
    #  import numpy
    #  import ctypes
    #  a1 = numpy.array((3, 1, 1, 1), dtype=numpy.int64)
    #  a2 = numpy.array((1, 1, 0, 0), dtype=numpy.int64)
    #  a1 & a2  # returns array([1, 1, 0, 0])
    #
    #  Load:
    #  i = numpy.array([1], dtype=numpy.int64)
    #  register = 0
    #  for refname in references:
    #      if refidx[refname] >= 63:
    #          register += 1
    #      a1[register] |= i[0] << (refidx[refname] - 63 * register)
    #  refmap[sample][chr] = (register, )
    #  
    # mask64 = ctypes.c_uint64(~0).value
    # 

    ############################################################################
    ## CLUSTERING PASS 1: EXACT-MATCH PATTERNS
    ############################################################################
    
    # distinct_pattern_index_map maps patterns to distinct_pattern_groups,
    # each containing indices in the ortho.groups table with that pattern
    distinct_pattern_index_map = dict()
    distinct_pattern_counts = list()
    distinct_pattern_groups = list()

    num_distinct_patterns = 0
    for i in range(num(ortho.groups)):
        pattern_incl_unanchored = tuple(map(tuple, ortho.groups[i]))
        pattern_excl_unanchored = [tuple()] * num_samples
        for sample in ortho.samples:
            pattern_excl_unanchored[sample.index] = \
                tuple(sorted(filter_unplaced_sequences(
                    pattern_incl_unanchored[sample.index],
                    config.samples[sample.id],
                    is_placed,
                    ignore_unplaced,
                    aggregate_unplaced=True
                )
            ))
        pattern = tuple(pattern_excl_unanchored)

        if pattern not in distinct_pattern_index_map:
            distinct_pattern_groups.append([])
            distinct_pattern_counts.append(0)
            distinct_pattern_index_map[pattern] = num_distinct_patterns
            num_distinct_patterns += 1

        distinct_pattern_groups[distinct_pattern_index_map[pattern]].append(i)
        distinct_pattern_counts[distinct_pattern_index_map[pattern]] += 1

        
    distinct_patterns_ranked = sorted(
        distinct_pattern_index_map,
        key=lambda p: distinct_pattern_counts[distinct_pattern_index_map[p]],
        reverse=True
    )

    distinct_patterns_excl_unanchored = [None] * num_distinct_patterns
    distinct_patterns_as_lists = [None] * num_distinct_patterns
    distinct_patterns_as_sets  = [None] * num_distinct_patterns
    distinct_pattern_clustered = [0] * num_distinct_patterns

    
    ############################################################################
    ## CLUSTERING PASS 2: MATCHING PATTERNS, ALLOWING MISSING DATA
    ############################################################################
    
    # clustered_pattern_index_map and partial_pattern_index_map map complete and
    # partial patterns, respectively, to clustered_pattern_groups and
    # partial_pattern_groups, each group containing indices of patterns in
    # distinct_patterns_ranked
    clustered_pattern_index_map = dict()
    clustered_pattern_representatives = list()
    clustered_pattern_groups = dict()
    clustered_pattern_counts = dict()

    partial_pattern_index_map = dict()
    partial_pattern_representatives = list()
    partial_pattern_counts = dict()
    partial_pattern_groups = dict()
    partial_patterns = []
    
    num_clustered_patterns = 0
    for i in range(num_distinct_patterns):
        distinct_patterns_excl_unanchored[i] = list(distinct_patterns_ranked[i])  # a tuples
        distinct_patterns_as_lists[i] = list(distinct_patterns_ranked[i])
        distinct_patterns_as_sets[i] = set()
        
        n = 0
        for j in range(num_samples):
            if distinct_patterns_as_lists[i][j]:
                distinct_patterns_excl_unanchored[i][j] = set(
                    distinct_patterns_excl_unanchored[i][j]
                )
                n += int(num(distinct_patterns_excl_unanchored[i][j]) > 0)
            else:
                distinct_patterns_as_lists[i][j] = set()
                distinct_patterns_excl_unanchored[i][j] = set()
                
            distinct_patterns_as_sets[i].update(distinct_patterns_as_lists[i][j])

        if n == num_samples:
            # patterns with complete chromosome membership are de facto their
            # own and complete clusters and require no further clustering.
            pattern = _as_tuples(distinct_patterns_excl_unanchored[i])
            if pattern not in clustered_pattern_index_map:
                clustered_pattern_index_map[pattern] = i  # sorted by rank
                clustered_pattern_groups[clustered_pattern_index_map[pattern]] = []
                clustered_pattern_counts[clustered_pattern_index_map[pattern]] = 0
                clustered_pattern_representatives.append(i)
            clustered_pattern_groups[clustered_pattern_index_map[pattern]].append(i)
            clustered_pattern_counts[clustered_pattern_index_map[pattern]] += \
                distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[i]]]
            distinct_pattern_clustered[i] = 1
        else:
            partial_patterns.append(i)

    clustered_pattern_representatives.sort(
        key=clustered_pattern_counts.get,
        reverse=True
    )

    print('Num input orthogroups:', num(ortho.groups), file=_STDERR)
    print('Num distinct patterns:', num_distinct_patterns, file=_STDERR)
    print('Num placed (samples complete):', num(clustered_pattern_representatives), file=_STDERR)
    print('Sum placed (samples complete):', sum(clustered_pattern_counts.values()), file=_STDERR)


    # Iterate through patterns at the bottom of the distinct_patterns_ranked
    # list and check if each pattern unambiguously matches a (more) complete
    # pattern (at the top of the list), and add it if it does.

    P = 0
    M = m = 0
    U = u = 0
    distinct_pattern_rep = dict()
    for j in range(num_distinct_patterns):  # a distinct_patterns_ranked index
        if distinct_pattern_clustered[j] != 0:
            continue
        
        best = list()
        for i in clustered_pattern_representatives:
            if num(distinct_patterns_as_sets[i].intersection(distinct_patterns_as_sets[j])) >= min_members:
                dist = calc_dist(
                    distinct_patterns_excl_unanchored[i],
                    distinct_patterns_excl_unanchored[j]
                )
                if dist < _EPSILON:
                    best.append(i)
                    
        if num(best) == 1:
            clustered_pattern_counts[best[0]] += distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]
            clustered_pattern_groups[best[0]].append(j)
            distinct_pattern_clustered[j] = 1
            distinct_pattern_rep[j] = best[0]
            P += 1

        elif num(best) > 1:
            # [M]ultiple matches
            m += 1
            M += distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]
            distinct_pattern_clustered[j] = -1
        else:
            # [U]nplaced
            u += 1
            U += distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]

            
    for j in range(num_distinct_patterns - 1, -1, -1):
        if distinct_pattern_clustered[j] != 0:
            continue

        best = list()
        for i in range(j):
            if num(distinct_patterns_as_sets[i].intersection(distinct_patterns_as_sets[j])) >= min_members:
                dist = calc_dist(
                    distinct_patterns_excl_unanchored[i],
                    distinct_patterns_excl_unanchored[j]
                )
                if dist < _EPSILON:
                    best.append(i)
                        
        if num(best) == 1:
            if distinct_pattern_clustered[best[0]] > 0:
                # pattern[i] was placed in another group
                while best[0] in distinct_pattern_rep:
                    best[0] = distinct_pattern_rep[best[0]]
            
            if best[0] not in clustered_pattern_groups:
                clustered_pattern_representatives.append(best[0])
                clustered_pattern_groups[best[0]] = [best[0]]
                clustered_pattern_counts[best[0]] = distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[best[0]]]]
                distinct_pattern_clustered[best[0]] = 1
            clustered_pattern_counts[best[0]] += distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]
            clustered_pattern_groups[best[0]].append(j)            
            distinct_pattern_clustered[j] = 1
            distinct_pattern_rep[j] = best[0]
            P += 1
        else:
            distinct_pattern_clustered[j] = -1
            

    ############################################################################
    ## CLUSTERING PASS 3: CHR-CHR ASSOCIATION CLUSTERING
    ############################################################################
    

    
    ############################################################################
    ## OUTPUT CLUSTERS
    ############################################################################
    
    print('Num placed (samples missing):', P, file=_STDERR)
    print('Sum placed (samples missing):', sum(partial_pattern_counts.values()), file=_STDERR)
    print('Num multiple best match:', m, file=_STDERR)
    print('Sum multiple best match:', M, file=_STDERR)
    print('Num unplaced (mismatches):', u, file=_STDERR)
    print('Sum unplaced (mismatches):', U, file=_STDERR)

    if cluster_counts_all_file or cluster_counts_mrg_file:
        if cluster_counts_mrg_file:
            cluster_counts_mrg_file.write(
                'Count\t%s\n' % OrthoFinderOrthogroups.format_orthogroups_header(ortho, id='Cluster')
            )
        if cluster_counts_all_file:
            cluster_counts_all_file.write(
                'Count\t%s\n' % OrthoFinderOrthogroups.format_orthogroups_header(ortho, id='Cluster')
            )
        cluster_count = 0
        for i in clustered_pattern_representatives:              # a distinct_patterns_ranked index
            for n, j in enumerate(clustered_pattern_groups[i]):  # a distinct_patterns_ranked index
                cluster_id = _CLUSTER_ID(cluster_count + 1)
                pattern_id = _PATTERN_ID(cluster_id, n + 1)
                
                num_patterns = num(distinct_pattern_groups[distinct_pattern_index_map[distinct_patterns_ranked[j]]])
                if num_patterns == 1:
                    pattern = ortho.groups[distinct_pattern_groups[distinct_pattern_index_map[distinct_patterns_ranked[j]]][0]]
                else:
                    pattern = distinct_patterns_ranked[j]
                
                if n == 0 and \
                   cluster_counts_mrg_file:
                    num_orthogroups = clustered_pattern_counts[i]                    
                    if min_orthogroups <= num_orthogroups <= max_orthogroups:
                        cluster_counts_mrg_file.write('%d\t%s\n' % (
                            num_orthogroups,
                            OrthoFinderOrthogroups.format_orthogroups_record(ortho, cluster_id, pattern)
                        ))
                        
                if cluster_counts_all_file:
                    num_orthogroups = distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]                
                    if min_orthogroups <= num_orthogroups <= max_orthogroups:
                        cluster_counts_all_file.write('%d\t%s\n' % (
                            num_orthogroups,
                            OrthoFinderOrthogroups.format_orthogroups_record(ortho, pattern_id, pattern)
                        ))
            cluster_count += 1
            
                        
        for j in range(num_distinct_patterns):
            if distinct_pattern_clustered[j] > 0:
                continue

            cluster_id = _CLUSTER_ID(cluster_count + 1)
            pattern_id = _PATTERN_ID(cluster_id, n + 1)
            
            num_orthogroups = distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]

            num_patterns = num(distinct_pattern_groups[distinct_pattern_index_map[distinct_patterns_ranked[j]]])
            if num_patterns == 1:
                pattern = ortho.groups[distinct_pattern_groups[distinct_pattern_index_map[distinct_patterns_ranked[j]]][0]]
            else:
                pattern = distinct_patterns_ranked[j]
            
            if min_orthogroups <= num_orthogroups <= max_orthogroups:
                if cluster_counts_mrg_file:
                    cluster_counts_mrg_file.write('%d\t%s\n' % (
                        distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]],
                        OrthoFinderOrthogroups.format_orthogroups_record(ortho, cluster_id, pattern)
                    ))
                if cluster_counts_all_file:
                    cluster_counts_all_file.write('%d\t%s\n' % (
                        distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]],
                        OrthoFinderOrthogroups.format_orthogroups_record(ortho, pattern_id, pattern)
                    ))            
                    
            cluster_count += 1
            
        if cluster_counts_mrg_file:
            cluster_counts_mrg_file.close()
        if cluster_counts_all_file:
            cluster_counts_all_file.close()

        
    if cluster_map_file:
        cluster_map_file.write('Cluster\t%s\n' % OrthoFinderOrthogroups.format_orthogroups_header(ortho))
        cluster_count = 0
        for i in clustered_pattern_representatives:
            for n, j in enumerate(clustered_pattern_groups[i]):
                cluster_id = _CLUSTER_ID(cluster_count + 1)
                pattern_id = _PATTERN_ID(cluster_id, n + 1)
                
                num_orthogroups = distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]
                
                if min_orthogroups <= num_orthogroups <= max_orthogroups:
                    for g in distinct_pattern_groups[distinct_pattern_index_map[distinct_patterns_ranked[j]]]:
                        cluster_map_file.write('%s\t%s\n' % (
                            pattern_id,
                            OrthoFinderOrthogroups.format_orthogroups_record(ortho, index=g)
                        ))
            cluster_count += 1

            
        for j in range(num_distinct_patterns):
            if distinct_pattern_clustered[j] > 0:
                continue

            cluster_id = _CLUSTER_ID(cluster_count + 1)
            pattern_id = _PATTERN_ID(cluster_id, n + 1)            
            
            num_orthogroups = distinct_pattern_counts[distinct_pattern_index_map[distinct_patterns_ranked[j]]]
            
            if min_orthogroups <= num_orthogroups <= max_orthogroups:
                for g in distinct_pattern_groups[distinct_pattern_index_map[distinct_patterns_ranked[j]]]:
                    cluster_map_file.write('%s\t%s\n' % (
                        pattern_id,
                        OrthoFinderOrthogroups.format_orthogroups_record(ortho, index=g)
                    ))
            cluster_count += 1
            
        cluster_map_file.close()
                                        
                
            
main(sys.argv[1:])



