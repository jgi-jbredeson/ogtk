
import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Join OrthoFinder Orthogroups.tsv files'

import sys
import getopt
import og.api.join as api

from og.api.join import ENUM_JOIN
from og.cli.utils import maxbits
from og.constants import _EMPTY, _EOL
from og.core.compression import STDIO
from og.core.orthogroups.parsers import OrthogroupsFile


num = len

                

def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    numbits = maxbits()
    
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <left.tsv> <right.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -f,--full-outer-join\n")
    stream.write("     Perform a full outer join [inner]\n")
    stream.write("\n")
    stream.write("  -I,--reverse-id-order\n")
    stream.write("     By default, the IDs of joined orthogroups are concatenated together\n")
    stream.write("     as `leftID:rightID`, this option instead outputs `rightID:leftID`\n")
    stream.write("\n")    
    stream.write("  -L,--max-intersecting-left <uint>\n")
    stream.write("     Maximum number of right orthogroups per left orthogroup [2^%d-1]\n" % numbits)
    stream.write("\n")    
    stream.write("  -l,--left-outer-join\n")
    stream.write("     Perform a left outer join [inner]\n")
    stream.write("\n")
    stream.write("  -m,--min-intersecting-members <uint>\n")
    stream.write("     Minimum number of intersecting members permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -M,--max-intersecting-members <uint>\n")
    stream.write("     Maximum number of intersecting members permitted per orthogroup [2^%d-1]\n" % numbits)
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -R,--max-intersecting-right <uint>\n")
    stream.write("     Maximum number of left orthogroup per right orthogroup [2^%d-1]\n" % numbits)
    stream.write("\n")
    stream.write("  -r,--right-outer-join\n")
    stream.write("     Perform a right outer join [inner]\n")
    stream.write("\n")
    stream.write("  -s,--min-intersecting-samples <uint>\n")
    stream.write("     Minimum number of intersecting samples permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -S,--max-intersecting-samples <uint>\n")
    stream.write("     Maximum number of intersecting samples permitted per orthogroup [2^%d-1]\n" % numbits)
    stream.write("\n")
    stream.write("  -x,--output-intersecting-strict\n")
    stream.write("     When performing the intersection, drop individual members from\n")
    stream.write("     left.tsv that are not also in right.tsv\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  1. left.tsv and right.tsv files are Orthogroups.tsv files with headers\n")
    stream.write("     defined. They must contain the same member types (e.g., loci or\n")
    stream.write("     references)\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)


    
def main(argv):
    short_flags = 'haAc:fL:lm:M:o:Q:R:rs:S:T:x'
    long_flags = (
        'help',
        'full-outer-join',
        'left-outer-join',
        'right-outer-join',
        'min-intersecting-members=',
        'max-intersecting-members=',
        'min-intersecting-samples=','min-intersecting-species=',
        'max-intersecting-samples=','max-intersecting-species=',
        'max-intersecting-left=',
        'max-intersecting-right=',
        'output-file=',
        'output-intersecting-members',
        'output-intersecting-strict',
        'reverse-id-order'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    method = ENUM_JOIN.INNER
    reverse_id_order = False
    output_file = STDIO
    output_intersecting_strict = False
    max_intersecting_left = sys.maxsize
    max_intersecting_right = sys.maxsize
    min_intersecting_members = 1
    max_intersecting_members = sys.maxsize
    min_intersecting_samples = 1
    max_intersecting_samples = sys.maxsize    
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-f','--full-outer-join'):
            method = ENUM_JOIN.FULL
        elif flag in ('-I','--reverse-id-order'):
            reverse_id_order = True            
        elif flag in ('-L','--max-left'):
            max_intersecting_left = int(value)
        elif flag in ('-l','--left-outer-join'):
            method = ENUM_JOIN.LEFT
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-m','--min-intersecting-members'):
            min_intersecting_members = int(value)
        elif flag in ('-M','--max-intersecting-members'):
            max_intersecting_members = int(value)
        elif flag in ('-R','--max-right'):
            max_intersecting_right = int(value)
        elif flag in ('-r','--right-outer-join'):
            method = ENUM_JOIN.RIGHT
        elif flag in ('-s',
                      '--min-intersecting-samples',
                      '--min-intersecting-species'):
            min_intersecting_samples = int(value)
        elif flag in ('-S',
                      '--max-intersecting-samples',
                      '--max-intersecting-species'):
            max_intersecting_samples = int(value)
        elif flag in ('-x','--output-intersecting-strict'):
            output_intersecting_strict = True

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    join(
        arguments[0], arguments[1], method=method,
        max_intersecting_left=max_intersecting_left,
        max_intersecting_right=max_intersecting_right,
        min_intersecting_members=min_intersecting_members,
        max_intersecting_members=max_intersecting_members,
        min_intersecting_samples=min_intersecting_samples,
        max_intersecting_samples=max_intersecting_samples,
        output_intersecting_strict=output_intersecting_strict,
        output_file=output_file,
        reverse_id_order=reverse_id_order,
    )

    return 0



def join(
        left_ortho_file,
        rght_ortho_file,
        method=ENUM_JOIN.INNER,
        max_intersecting_left=sys.maxsize,
        max_intersecting_right=sys.maxsize,
        min_intersecting_members=1,
        max_intersecting_members=sys.maxsize,
        min_intersecting_samples=1,
        max_intersecting_samples=sys.maxsize,
        output_intersecting_strict=False,
        output_file=sys.stdout,
        reverse_id_order=False
):
    left_ortho = OrthogroupsFile(left_ortho_file)
    rght_ortho = OrthogroupsFile(rght_ortho_file)        

    join_ortho = \
        api.join(
            left_ortho, rght_ortho, method=method,
            max_intersecting_left=max_intersecting_left,
            max_intersecting_right=max_intersecting_right,
            min_intersecting_members=min_intersecting_members,
            max_intersecting_members=max_intersecting_members,
            min_intersecting_samples=min_intersecting_samples,
            max_intersecting_samples=max_intersecting_samples,
            output_intersecting_strict=output_intersecting_strict,
            reverse_id_order=reverse_id_order,
        )

    join_ortho.to_file(output_file)
