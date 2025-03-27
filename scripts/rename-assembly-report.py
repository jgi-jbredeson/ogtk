#!/usr/bin/env python3

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Rename assembly_report sequences'

import re
import sys
import string
import getopt

from og.core.parsers.assembly_report import AssemblyReportFile
from og.core.parsers.assembly_report import (
    _VALID_ATTRIBUTES,
    _VALID_COLUMNS,
    UNIT_PRIMARY,
)
from og.constants import (
    _COLON,
    _EMPTY,
    _TAB,
    dict
)

_ALPHA = str.maketrans({a:None for a in string.ascii_letters})

_RNAME_COLUMNS = ('Sequence-Name','GenBank-Accn','RefSeq-Accn','UCSC-style-name')
_COLUMN_ATTR_MAP = dict(zip(_VALID_COLUMNS, _VALID_ATTRIBUTES))

num = len


class Pattern(object):
    def __init__(self, column=None, pattern=None):
        self.column = column
        self.pattern = pattern

        
def get_pattern(value):
    fields = value.split(_COLON,maxsplit=1)
    if num(fields) == 2:
        if fields[0] in _VALID_COLUMNS:
            if fields[0] in _RNAME_COLUMNS:
                return Pattern(_COLUMN_ATTR_MAP[fields[0]], re.compile(fields[1]))
            else:
                usage('Cannot write to column: not a name/accession field: ' + fields[0])
        else:
            usage('Invalid column name (case sensitive): ' + fields[0])
    else:
        return Pattern(_COLUMN_ATTR_MAP[_RNAME_COLUMNS[0]], re.compile(value))


def get_subgroup(pattern, field):
    match = pattern.search(field)
    if match:
        for group in match.groups():
            if group is not None:
                return group
    else:
         usage('Unmatched pattern: ' + field)
    return field


def format_number(field, width=0):
    number = field.translate(_ALPHA)
    number_width = len(number)
    target_width = 0 if width < 0 else width
    if number.isdigit() and \
       number_width < target_width:
        return ('0' * (target_width - number_width)) + field
    else:
        return field

    
def _format_id(fields, column_indices, count_dict, regexp=None, zwidth=0):
    if regexp is None:
        if fields[2] in count_dict:
            count_dict[fields[2]] += 1
        else:
            count_dict[fields[2]] = 1
        label = str(count_dict[fields[2]])
    else:
        label = get_subgroup(
            regexp.pattern,
            fields[column_indices[regexp.column]]
        )
    return format_number(label, zwidth)


def format_id(record, count_dict, regexp=None, zwidth=0):
    if regexp is None:
        if record.assigned_molecule in count_dict:
            count_dict[record.assigned_molecule] += 1
        else:
            count_dict[record.assigned_molecule] = 1
        label = str(count_dict[record.assigned_molecule])
    else:
        label = get_subgroup(
            regexp.pattern,
            getattr(record, regexp.column)
        )
    return format_number(label, zwidth)
    
    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <assembly_report.txt> <prefix>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0         10         20       30        40        50        60        70        80
    stream.write("  -a,--alt-scaffold-regexp [column-name:]<pattern>\n")
    stream.write("     Use specified regular expression (which must contain at least one\n")
    stream.write("     capture group) to parse an unique identifier from the value given\n")
    stream.write("     in column-name for alt-scaffold sequences.\n")
    stream.write("\n")
    stream.write("  -f,--fix-patch-regexp [column-name:]<pattern>\n")
    stream.write("     Use specified regular expression (which must contain at least one\n")
    stream.write("     capture group) to parse an unique identifier from the value given\n")
    stream.write("     in column-name for fix-patch sequences.\n")
    stream.write("\n")
    stream.write("  -n,--novel-patch-regexp [column-name:]<pattern>\n")
    stream.write("     Use specified regular expression (which must contain at least one\n")
    stream.write("     capture group) to parse an unique identifier from the value given\n")
    stream.write("     in column-name for novel-patch sequences.\n")
    stream.write("\n")
    stream.write("  -U,--unlocalized-scaffold-regexp [column-name:]<pattern>\n")
    stream.write("     Use specified regular expression (which must contain at least one\n")
    stream.write("     capture group) to parse an unique identifier from the value given\n")
    stream.write("     in column-name for unlocalized-scaffold sequences.\n")
    stream.write("\n")    
    stream.write("  -u,--unplaced-scaffold-regexp [column-name:]<pattern>\n")
    stream.write("     Use specified regular expression (which must contain at least one\n")
    stream.write("     capture group) to parse an unique identifier from the value given\n")
    stream.write("     in column-name for unplaced-scaffold sequences.\n")
    stream.write("\n")    
    stream.write("  -Z,--chromosome-zero-pad-width <uint>\n")
    stream.write("     Generate zero-padded chromosome numbers with the specified width.\n")
    stream.write("     Not applied to lettered sex chromosomes.\n")
    stream.write("     (default: 0)\n")
    stream.write("\n")
    stream.write("  -z,--scaffold-zero-pad-width <uint>\n")
    stream.write("     Generate zero-padded scaffold numbers with the specified width.\n")
    stream.write("     Also applies to the values of the above regexp captured groups.\n")
    stream.write("     (default: 0)\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit.\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0         10         20       30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  - By default, the sequence names for all sequence-roles, except\n")
    stream.write("    assembled-molecule, are re-numbered.\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)
    

    

def main(argv):
    short_flags = 'Z:z:U:u:a:n:f:h'
    long_flags = (
        'chromosome-zero-pad-width=',
        'scaffold-zero-pad-width=',
        'unlocalized-scaffold-regexp=',
        'unlaced-scaffold-regexp=',
        'alt-scaffold-regexp=',
        'novel-patch-regexp=',
        'fix-patch-regexp=',
        'help'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    chr_zero_pad_width = 1
    sca_zero_pad_width = 1        
    fix_patch_regexp = None
    novel_patch_regexp = None
    alt_scaffold_regexp = None
    unplaced_scaffold_regexp = None
    unlocalized_scaffold_regexp = None
    for flag, value in options:
        if flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-a','--alt-scaffold-regexp'):
            alt_scaffold_regexp = get_pattern(value)
        elif flag in ('-U','--unlocalized-scaffold-regexp'):
            unlocalized_scaffold_regexp = get_pattern(value)
        elif flag in ('-u','--unplaced-scaffold-regexp'):
            unplaced_scaffold_regexp = get_pattern(value)
        elif flag in ('-n','--novel-patch-regexp'):
            novel_patch_regexp = get_pattern(value)
        elif flag in ('-f','--fix-patch-regexp'):
            fix_patch_regexp = get_pattern(value)
        elif flag in ('-Z','--chromosome-zero-pad-width'):
            chr_zero_pad_width = int(value)
        elif flag in ('-z','--scaffold-zero-pad-width'):
            sca_zero_pad_width = int(value)            

    if num(arguments) != 2:
        usage('Unexpected number of arguments')
        
    assembly_report = AssemblyReportFile(arguments[0])
    assembly_prefix = arguments[1]

    alt_count = dict()
    fix_count = dict()
    novel_count = dict()
    unplaced_count = 0
    unlocalized_count = dict()
    print('# ' + _TAB.join(_VALID_COLUMNS))
    
    for record in assembly_report.values():
        if record.is_primary:
            if record.assigned_type == 'Chromosome':
                if record.is_assembled_molecule:
                    record.sequence_name \
                        = assembly_prefix \
                        + format_number(
                            record.assigned_molecule,
                            chr_zero_pad_width
                        )

                elif record.is_unlocalized:
                    record.sequence_name = '%s%s.Un%s' % (
                        assembly_prefix,
                        format_number(
                            record.assigned_molecule,
                            chr_zero_pad_width
                        ),
                        format_id(
                            record,
                            unlocalized_count,
                            unlocalized_scaffold_regexp,
                            sca_zero_pad_width
                        )
                    )
                else:
                    raise NotImplementedError('%s:%s:%s' % (
                        record.assigned_type,
                        record.sequence_role,
                        record.sequence_name
                    ))
            
            elif record.is_unplaced:
                if unplaced_scaffold_regexp:
                    record.sequence_name = get_subgroup(
                        unplaced_scaffold_regexp.pattern,
                        getattr(record, unplaced_scaffold_regexp.column)
                    )
                    record.sequence_name = '%sUn%s' % (
                        assembly_prefix,
                        format_number(
                            record.sequence_name,
                            sca_zero_pad_width
                        )
                    )
                else:
                    unplaced_count += 1
                    record.sequence_name = '%sUn%s' % (
                        assembly_prefix,
                        format_number(
                            str(unplaced_count),
                            sca_zero_pad_width
                        )
                    )
            else:
                raise NotImplementedError('%s:%s:%s' % (
                    record.assigned_type,
                    record.sequence_role,
                    record.sequence_name
                ))


        elif record.is_nonnuclear:
            if record.is_assembled_molecule:
                molecule = record.assigned_type
                if molecule == 'chloroplast':
                    # Assigned-Molecule is 'Pltd', not specific enough.
                    record.sequence_name = assembly_prefix + 'CP'
                else:
                    # Plastid
                    # Mitochondrion
                    # Mitochondrial Plasmid
                    record.sequence_name = assembly_prefix + record.assigned_molecule
            else:
                raise NotImplementedError('%s:%s:%s' % (
                    record.assigned_type,
                    record.sequence_role,
                    record.sequence_name
                ))
            
        else:  # Assembly-Unit not 'Primary Assembly'
            if record.assigned_type == 'Chromosome':
                if record.is_alternate:
                    record.sequence_name = '%s%s.Alt%s' % (
                        assembly_prefix,
                        format_number(
                            record.assigned_molecule,
                            chr_zero_pad_width
                        ),
                        format_id(
                            record,
                            alt_count,
                            alt_scaffold_regexp,
                            sca_zero_pad_width
                        )
                    )
                elif record.is_patch:
                    patch_type = \
                        'NP' if record.sequence_role == 'novel-patch' else \
                        'FP' if record.sequence_role == 'fix-patch' else \
                        'OP' # other-patch?
                        
                    record.sequence_name = '%s%s.%s%s' % (
                        assembly_prefix,
                        format_number(
                            record.assigned_molecule,
                            chr_zero_pad_width
                        ),
                        patch_type,
                        format_id(
                            record,
                            novel_count,
                            novel_patch_regexp,
                            sca_zero_pad_width
                        )
                    )
                else:
                    raise NotImplementedError('%s:%s:%s' % (
                        record.assigned_type,
                        record.sequence_role,
                        record.sequence_name
                    ))
            else:
                raise NotImplementedError('%s:%s:%s' % (
                    record.assigned_type,
                    record.sequence_role,
                    record.sequence_name
                ))

        print(str(record))


if __name__ == '__main__':
    main(sys.argv[1:])

    
