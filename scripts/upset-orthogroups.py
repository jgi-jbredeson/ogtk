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
from og.core.members import map_loci_to_sequence_counts
from og.core.members import filter_unplaced_sequences
from og.core.compression import open, STDIO
from og.core.parsers.config import SampleConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
from og.constants import _TAB, _SPACE, _EMPTY, _EOL


num = len
_GROUP_MEMBERSHIP = 0x1
_GROUP_MULTIPLES = 0x2
_SORT_CLUSTERS = 0x1
_SORT_MEMBERS = 0x2
_MULTI = '\u25CF'
_MINUS = '\u2015'

_UNICODE_NUMMAP = {
    '1' :'\u2460',
    '2' :'\u2461',
    '3' :'\u2462',
    '4' :'\u2463',
    '5' :'\u2464',
    '6' :'\u2465',
    '7' :'\u2466',
    '8' :'\u2467',
    '9' :'\u2468',
    # '10':'\u2469',
    # '11':'\u246A',
    # '12':'\u246B',
    # '13':'\u246C',
    # '14':'\u246D',
    # '15':'\u246E',
    # '16':'\u246F',
    # '17':'\u2470',
    # '18':'\u2471',
    # '19':'\u2472',
    # '20':'\u2473',
}


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
    stream.write("  -a,--force-ascii\n")
    stream.write("     Force writing output in ASCII-only characters\n")
    stream.write("\n")
    stream.write("  -g,--group-by <enum>\n")
    stream.write("     Group output by (1) membership, (2) number of multiples\n")
    stream.write("     (See the `--max-count` option below), or perform no grouping [0]\n")
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
    stream.write("     Count up to `-M` number of sequences per sample and orthogroup.\n")
    stream.write("     Counts greater than this threshold are converted to the `%s` character\n" % _MULTI)
    stream.write("     (or `*` if `--force-ascii` is enabled).\n")
    stream.write("\n")
    stream.write("  -n,--map-to-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then calculate plot.\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -p,--input-sequence-names\n")
    stream.write("     Locus names have already been mapped to their corresponding sequence\n")
    stream.write("     names in the input orthogroups file. Perform filtering accordingly.\n")
    stream.write("\n")
    stream.write("  -s,--sort-by <enum>\n")
    stream.write("     Sort plot rows by (1) number of clusters or (2) number of members [-1]\n")
    stream.write("     (Negate the enumerative value to reverse sort order)\n")
    stream.write("\n")
    stream.write("  -S,--upset-separator <str>\n")
    stream.write("     Add space between columns of the upset plot using the given separator\n")
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
    short_flags = 'ahS:s:g:M:o:pniIu'
    long_flags = (
        'force-ascii',
        'help',
        'upset-sep=',
        'sort-by=',
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

    sort_by = -1
    group_by = 0
    max_count = 0
    is_placed = _placed
    map_seq_names = False
    input_seq_names = False
    output_file = STDIO
    ignore_unplaced = 0
    sep = _EMPTY
    multi = _MULTI
    minus = _MINUS
    force_ascii = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-a','--force-ascii'):
            force_ascii = True
            multi = '*'
            minus = '-'
            sep = _SPACE
        elif flag in ('-g','--group-by'):
            group_by |= (0x1 | int(value))
        elif flag in ('-s','--sort-by'):
            sort_by = int(value)
        elif flag in ('-S','--upset-separator'):
            sep = value.encode('utf-8').decode('unicode_escape')
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
            is_placed = _localized

    if not (abs(sort_by) & (_SORT_CLUSTERS | _SORT_MEMBERS)):
        usage("Invalid enumerative value: --sort-by=%s" % str(sort_by))
        
    if not (input_seq_names or map_seq_names):
        output_seq_names = False
        ignore_unplaced = False            
            
    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    ortho  = OrthoFinderOrthogroups(arguments[0])
    config = SampleConfigFile(arguments[1], load_files=True, map_assigned_molecule=True)
    output = open(output_file, 'w')
    
    for sample in ortho.samples:
        if sample.id not in config.samples:
            raise KeyError(
                "Sample not found in YAML file: '%s'" % str(sample.id)
            )

    member_counts = dict()
    pattern_counts = dict()
    num_samples = num(ortho.samples)
    for group in ortho.groups:
        num_members = 0
        pattern = [0] * num_samples
        sequence_names = [None] * num_samples
        if map_seq_names:
            for sample in ortho.samples:
                sequence_names[sample.index] = \
                    map_loci_to_sequence_counts(
                        group[sample.index],
                        config.samples[sample.id],
                        is_placed,
                        ignore_unplaced=False
                    )
        else:  # either loci or pre-mapped sequences:
            for sample in ortho.samples:
                sequence_names[sample.index] = \
                    count_items(
                        group[sample.index]
                    )
                 

        for sample in ortho.samples:
            if map_seq_names or input_seq_names:
                members = filter_unplaced_sequences(
                    sequence_names[sample.index],
                    config.samples[sample.id],
                    is_placed,
                    ignore_unplaced,
                    aggregate_unplaced=False
                )
            else:
                members = group[sample.index]

            n = num(members)
            
            pattern[sample.index] = minus if n < 1 else multi if n > max_count else str(n)
            num_members += n

        try:
            member_counts[tuple(pattern)] += num_members
            pattern_counts[tuple(pattern)] += 1
        except KeyError:
            member_counts[tuple(pattern)] = num_members
            pattern_counts[tuple(pattern)] = 1
            
    output.write(
        _TAB.join((
            _SPACE.join(s.id for s in ortho.samples),
            'Clusters',
            'References' if (input_seq_names or map_seq_names) else 'Proteins'
        )) + _EOL
    )

    if force_ascii or max_count > num(_UNICODE_NUMMAP):
        sep = _SPACE if sep is _EMPTY else sep
        def _fmtchar(x):
            return str(x)
    else:
        def _fmtchar(x):
            try:
                return _UNICODE_NUMMAP[x]
            except KeyError:
                return str(x)

    if abs(sort_by) & _SORT_CLUSTERS:
        def _sort(x):
            return (pattern_counts.get(x, 0), member_counts.get(x,0))
    elif abs(sort_by) & _SORT_MEMBERS:
        def _sort(x):
            return (member_counts.get(x,0), pattern_counts.get(x, 0))
    else:
        raise AssertionError('Invalid sort value: %s' % str(sort_by))
        
        
    if group_by & _GROUP_MEMBERSHIP:
        group_patterns = dict()
        group_counts = dict()
        if group_by & _GROUP_MULTIPLES:
            for pattern in pattern_counts:
                _pattern = tuple(sorted(pattern))
                _nmulti = _pattern.count(multi)
                _counts = _sort(pattern)
                if _nmulti not in group_patterns:
                    group_patterns[_nmulti] = dict()
                    group_counts[_nmulti] = [0,0]
                group_patterns[_nmulti][pattern] = _counts
                group_counts[_nmulti][0] += _counts[0]
                group_counts[_nmulti][1] += _counts[1]
        else:
            for pattern in pattern_counts:
                _pattern = tuple(sorted(pattern))
                _counts = _sort(pattern)
                if _pattern not in group_patterns:
                    group_patterns[_pattern] = dict()
                    group_counts[_pattern] = [0,0]
                group_patterns[_pattern][pattern] = _counts
                group_counts[_pattern][0] += _counts[0]
                group_counts[_pattern][1] += _counts[1]

        for _pattern in sorted(group_patterns, key=group_counts.get, reverse=(sort_by < 0)):
            output.write("##total=%d\n" % group_counts[_pattern][0])
            for pattern in sorted(group_patterns[_pattern], key=group_patterns[_pattern].get, reverse=(sort_by < 0)):
                output.write(
                    _TAB.join((
                        sep.join(map(_fmtchar, pattern)),
                        str(pattern_counts[pattern]),
                        str(member_counts[pattern])
                    )) + _EOL
                )
    else:
        for pattern in sorted(pattern_counts, key=_sort, reverse=(sort_by < 0)):
            output.write(
                _TAB.join((
                    sep.join(map(_fmtchar, pattern)),
                    str(pattern_counts[pattern]),
                    str(member_counts[pattern])
                )) + _EOL
            )
            
    output.close()

    return 0
    

if __name__ == '__main__':
    main(sys.argv[1:])

