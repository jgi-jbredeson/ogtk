
from og import (
    __authors__,
    __contact__,
    __pkgname__,
    __version__,
    __program__,
)
__purpose__ = 'Subset OrthoFinder Orthogroups.tsv file members'

import os
import sys
import getopt
import og.api.subset as api

from og.constants import _EMPTY
from og.cli.utils import read_list, get_exe
from og.api.subset import ENUM_CALLKEY
from og.core.compressio import STDIO
from og.core.orthogroups.parsers import OrthogroupsFile



num = len



def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <in.tsv>\n" % get_exe(__program__, __file__))
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -g,--include-groups <list>\n")
    stream.write("     Subset orthogroups by ID, keeping only the groups listed\n")
    stream.write("\n")
    stream.write("  -G,--exclude-groups <list>\n")
    stream.write("     Subset orthogroups by ID, keeping only the groups not listed\n")
    stream.write("\n")
    stream.write("  -m,--include-members <list>\n")
    stream.write("     Subset orthogroup members by ID, keeping only the members listed\n")
    stream.write("\n")
    stream.write("  -M,--exclude-members <list>\n")
    stream.write("     Subset orthogroup members by ID, keeping only the members not listed\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -s,--include-samples <list>\n")
    stream.write("     Subset orthogroup samples by ID, keeping only samples listed\n")
    stream.write("\n")
    stream.write("  -S,--exclude-samples <list>\n")
    stream.write("     Subset orthogroup samples by ID, keeping only samples not listed\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    stream.write("\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("Notes:\n")
    stream.write("  1. The in.tsv file is an Orthogroups.tsv file with header defined.\n")
    stream.write("\n")
    stream.write("  2. Input list arguments may be comma-separated strings of IDs or the name\n")
    stream.write("     of a file listing IDs, one per line.\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    
    
def main(argv):
    short_flags = 'ho:g:G:m:M:s:S:'
    long_flags = (
        'help',
        'output-file=',
        'include-groups=',
        'exclude-groups=',
        'include-members=',
        'exclude-members=',
        'include-samples=',
        'exclude-samples='
    )

    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    call_order = []
    output_file = STDIO
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-g','--include-groups'):
            call_order.append((
                ENUM_CALLKEY.INCLUDE_GROUPS,
                set(read_list(value))
            ))
        elif flag in ('-G','--exclude-groups'):
            call_order.append((
                ENUM_CALLKEY.EXCLUDE_GROUPS,
                set(read_list(value))
            ))
        elif flag in ('-m','--include-members'):
            call_order.append((
                ENUM_CALLKEY.INCLUDE_MEMBERS,
                set(read_list(value))
            ))
        elif flag in ('-M','--exclude-members'):
            call_order.append((
                ENUM_CALLKEY.EXCLUDE_MEMBERS,
                set(read_list(value))
            ))
        elif flag in ('-s','--include-samples'):
            call_order.append((
                ENUM_CALLKEY.INCLUDE_SAMPLES,
                set(read_list(value))
            ))
        elif flag in ('-S','--exclude-samples'):
            call_order.append((
                ENUM_CALLKEY.EXCLUDE_SAMPLES,
                set(read_list(value))
            ))
        
    if num(arguments) == 0:
        usage()
    if num(arguments) != 1:
        usage('Unexpected number of arguments')

    subset(arguments[0], call_order, output_file=output_file)

    return 0



def subset(ortho_file, call_order, output_file=sys.stdout):
    ortho = api.subset(OrthogroupsFile(ortho_file), call_order)
    ortho.to_file(output_file)
