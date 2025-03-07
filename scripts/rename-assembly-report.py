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

_ALPHA = str.maketrans({a:None for a in string.ascii_letters})

num = len

def _istype(value, typ):
    try:
        typ(value)
    except ValueError:
        return False
    return True



def format_number(field, width=0):
    number = field.translate(_ALPHA)
    number_width = len(number)
    target_width = 0 if width < 0 else width
    if number.isdigit() and \
       number_width < target_width:
        return ('0' * (target_width - number_width)) + field
    else:
        return field


    
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
    stream.write("  -a,--alt-scaffold-regexp <pattern>\n")
    stream.write("     Use this regular expression to remove leading characters from the\n")
    stream.write("     names of alt-scaffold sequences not associated with chromosomes.\n")
    stream.write("\n")
    stream.write("  -n,--novel-patch-regexp <pattern>\n")
    stream.write("     Use this regular expression to remove leading characters from the\n")
    stream.write("     names of novel-patch sequences not associated with chromosomes.\n")
    stream.write("\n")
    stream.write("  -f,--fix-patch-regexp <pattern>\n")
    stream.write("     Use this regular expression to remove leading characters from the\n")
    stream.write("     names of fix-patch sequences not associated with chromosomes.\n")
    stream.write("\n")
    stream.write("  -Z,--chromosome-zero-pad-width <uint>\n")
    stream.write("     Generate zero-padded chromosome numbers with the specified width.\n")
    stream.write("     Does not apply to lettered sex chromosomes.\n")
    stream.write("     (default: 0)\n")
    stream.write("\n")
    stream.write("  -z,--scaffold-zero-pad-width <uint>\n")
    stream.write("     Generate zero-padded scaffold numbers with the specified width\n")
    stream.write("     (default: 0)\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit.\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0         10         20       30        40        50        60        70        80
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)
    

    

def main(argv):
    short_flags = 'Z:z:a:n:f:h'
    long_flags = (
        'chromosome-zero-pad-width=',
        'scaffold-zero-pad-width=',
        'alt-scaffold-regexp=',
        'novel-patch-regexp=',
        'fix-patch-regexp=',
        'help'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)


    fix_patch_regexp = None
    novel_patch_regexp = None
    alt_scaffold_regexp = None
    chr_zero_pad_width = 1
    sca_zero_pad_width = 1
    for flag, value in options:
        if flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-a','--alt-scaffold-regexp'):
            alt_scaffold_regexp = re.compile(value)
        elif flag in ('-n','--novel-patch-regexp'):
            novel_patch_regexp = re.compile(value)
        elif flag in ('-f','--fix-patch-regexp'):
            fix_patch_regexp = re.compile(value)
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
    for line in assembly_report:
        line   = line.strip()
        fields = line.split('\t')
        
        if line.startswith('#'):
            print(line)
            continue

        if fields[7].lower() == 'primary assembly':   # Chromosome
            if fields[3].lower() == 'chromosome':            
                if fields[1] == 'assembled-molecule':
                    fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)
                
                elif fields[1] == 'unlocalized-scaffold':
                    fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)
                    
                    if fields[2] in unlocalized_count:
                        unlocalized_count[fields[2]] += 1
                    else:
                        unlocalized_count[fields[2]] = 1

                    fields[0] = '%s.Un%s' % (
                        fields[0],
                        format_number(str(unlocalized_count[fields[2]]), sca_zero_pad_width),
                    )
                else:
                    raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0]))
                
            elif fields[1] == 'unplaced-scaffold':
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
                    fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)

                    if fields[2] in alt_count:
                        alt_count[fields[2]] += 1
                    else:
                        alt_count[fields[2]] = 1

                    fields[0] += '.Alt%s' % (
                        format_number(str(alt_count[fields[2]]), sca_zero_pad_width),
                    )
                elif fields[1] == 'novel-patch':
                    fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)

                    if fields[2] in novel_count:
                        novel_count[fields[2]] += 1
                    else:
                        novel_count[fields[2]] = 1

                    fields[0] += '.NP%s' % (
                        format_number(str(novel_count[fields[2]]), sca_zero_pad_width),
                    )
                elif fields[1] == 'fix-patch':
                    fields[0] = assembly_prefix + format_number(fields[2], chr_zero_pad_width)

                    if fields[2] in fix_count:
                        fix_count[fields[2]] += 1
                    else:
                        fix_count[fields[2]] = 1

                    fields[0] += '.FP%s' % (
                        format_number(str(fix_count[fields[2]]), sca_zero_pad_width),
                    )
                else:
                    raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0]))
                    
            elif fields[1] == 'alt-scaffold':
                if alt_scaffold_regexp:
                    fields[0] = alt_scaffold_regexp.sub(r'',fields[0])
                fields[0] = assembly_prefix + fields[0]

            elif fields[1] == 'novel-patch':
                if novel_patch_regexp:
                    fields[0] = novel_patch_regexp.sub(r'',fields[0])
                fields[0] = assembly_prefix + fields[0]
                
            elif fields[1] == 'fix-patch':
                if fix_patch_regexp:
                    fields[0] = fix_patch_regexp.sub(r'',fields[0])
                fields[0] = assembly_prefix + fields[0]
                
            else:
                raise NotImplementedError('%s:%s:%s' % (fields[7],fields[1],fields[0])) 

        print('\t'.join(fields))


if __name__ == '__main__':
    main(sys.argv[1:])

    
