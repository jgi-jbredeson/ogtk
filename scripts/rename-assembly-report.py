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

from og.core.utils import index_list

_ALPHA = str.maketrans({a:None for a in string.ascii_letters})
_VALID_COLUMNS = ('Sequence-Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank-Accn','Relationship','RefSeq-Accn','Assembly-Unit','Sequence-Length','UCSC-style-name')
_RNAME_COLUMNS = ('Sequence-Name','GenBank-Accn','RefSeq-Accn','UCSC-style-name')

num = len


class Pattern(object):
    def __init__(self, column=None, pattern=None):
        self.column = column
        self.pattern = pattern

        
def get_pattern(value):
    fields = value.split(':',maxsplit=1)
    if num(fields) == 2:
        if fields[0] in _VALID_COLUMNS:
            if fields[0] in _RNAME_COLUMNS:
                return Pattern(fields[0], re.compile(fields[1]))
            else:
                usage('Cannot write to column: not a name/accession field: ' + fields[0])
        else:
            usage('Invalid column name: ' + fields[0])
    else:
        return Pattern(_RNAME_COLUMNS[0], re.compile(value))


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

    
def format_id(fields, column_indices, count_dict, regexp=None, zwidth=0):
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
            
    
    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = '' if message is None else 'ERROR: %s\n\n' % message
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
        
            
    assembly_report = open(arguments[0], 'r')
    assembly_prefix = arguments[1]

    alt_count = dict()
    fix_count = dict()
    novel_count = dict()
    unplaced_count = 0
    unlocalized_count = dict()
    column_indices = None
    for line in assembly_report:
        line   = line.strip()
        fields = line.split('\t')
        
        if line.startswith('#'):
            if line.startswith('# %s' % _VALID_COLUMNS[0]):
                fields[0] = fields[0].lstrip('#').strip()
                column_indices = index_list(fields)
            print(line)
            continue
        elif column_indices is None:
            print('# %s' % '\t'.join(_VALID_COLUMNS))
            column_indices = index_list(_VALID_COLUMNS)

        if fields[7].lower() == 'primary assembly':
            if fields[3].lower() == 'chromosome':
                if fields[1] == 'assembled-molecule':
                    fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)
                
                elif fields[1] == 'unlocalized-scaffold':
                    fields[0] = '%s%s.Un%s' % (
                        assembly_prefix,
                        format_number(fields[2], chr_zero_pad_width),
                        format_id(
                            fields,
                            column_indices,
                            unlocalized_count,
                            unlocalized_scaffold_regexp,
                            sca_zero_pad_width
                        )
                    )
                else:
                    raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0]))
                
            elif fields[1] == 'unplaced-scaffold':
                if unplaced_scaffold_regexp:
                    fields[0] = get_subgroup(
                        unplaced_scaffold_regexp.pattern,
                        fields[column_indices[unplaced_scaffold_regexp.column]]
                    )
                    fields[0] = '%sUn%s' % (
                        assembly_prefix,
                        format_number(fields[0], sca_zero_pad_width)
                    )
                else:
                    unplaced_count += 1
                    fields[0] = '%sUn%s' % (
                        assembly_prefix,
                        format_number(str(unplaced_count), sca_zero_pad_width)
                    )

            else:
                raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0]))

        elif fields[7] == 'non-nuclear':
            if fields[1] == 'assembled-molecule':
                molecule = fields[3].lower()
                if molecule == 'chloroplast':
                    # field[2] is 'Pltd', not specific enough.
                    fields[0] = assembly_prefix + 'CP'
                else:
                    # Plastid
                    # Mitochondrion
                    # Mitochondrial Plasmid
                    fields[0] = assembly_prefix + fields[2]
            else:
                raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0]))               
            
        else:  # not 'Primary Assembly'
            if fields[3].lower() == 'chromosome':            
                if fields[1] == 'alt-scaffold':
                    fields[0] = '%s%s.Alt%s' % (
                        assembly_prefix,
                        format_number(fields[2], chr_zero_pad_width),
                        format_id(
                            fields,
                            column_indices,
                            alt_count,
                            alt_scaffold_regexp,
                            sca_zero_pad_width
                        )
                    )
                    # fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)

                    # if alt_scaffold_regexp:
                    #     label = get_subgroup(
                    #         alt_scaffold_regexp.pattern,
                    #         fields[column_indices[alt_scaffold_regexp.column]]
                    #     )
                    # else:
                    #     if fields[2] in alt_count:
                    #         alt_count[fields[2]] += 1
                    #     else:
                    #         alt_count[fields[2]] = 1
                    #     label = str(alt_count[fields[2]])

                    # fields[0] += '.Alt%s' % (
                    #     format_number(label, sca_zero_pad_width),
                    # )
                elif fields[1] == 'novel-patch':
                    fields[0] = '%s%s.NP%s' % (
                        assembly_prefix,
                        format_number(fields[2], chr_zero_pad_width),
                        format_id(
                            fields,
                            column_indices,
                            novel_count,
                            novel_patch_regexp,
                            sca_zero_pad_width
                        )
                    )
                    # fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)

                    # if novel_patch_regexp:
                    #     label = get_subgroup(
                    #         novel_patch_regexp.pattern,
                    #         fields[column_indices[novel_patch_regexp.column]]
                    #     )
                    # else:
                    #     if fields[2] in novel_count:
                    #         novel_count[fields[2]] += 1
                    #     else:
                    #         novel_count[fields[2]] = 1
                    #     label = str(novel_count[fields[2]])
                        
                    # fields[0] += '.NP%s' % (
                    #     format_number(label, sca_zero_pad_width),
                    # )
                elif fields[1] == 'fix-patch':
                    fields[0] = '%s%s.FP%s' % (
                        assembly_prefix,
                        format_number(fields[2], chr_zero_pad_width),
                        format_id(
                            fields,
                            column_indices,
                            fix_count,
                            fix_patch_regexp,
                            sca_zero_pad_width
                        )
                    )
                    
                    # fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)

                    # if fix_patch_regexp:
                    #     label = get_subgroup(
                    #         fix_patch_regexp.pattern,
                    #         fields[column_indices[fix_patch_regexp.column]]
                    #     )
                    # else:
                    #     if fields[2] in fix_count:
                    #         fix_count[fields[2]] += 1
                    #     else:
                    #         fix_count[fields[2]] = 1
                    #     label = str(fix_count[fields[2]])
                        
                    # fields[0] += '.FP%s' % (
                    #     format_number(label, sca_zero_pad_width),
                    # )
                else:
                    raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0]))
            else:
                raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0])) 

        print('\t'.join(fields))


if __name__ == '__main__':
    main(sys.argv[1:])

    
