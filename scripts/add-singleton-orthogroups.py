
import os
import sys
import getopt

from math import inf as _POS_INF
from og.core.utils import index_list
from og.core.members import map_loci_to_sequence_counts
from og.core.parsers.config import SampleConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
from og.constants import (
    _COLON,
    _COMMENT,
    _EMPTY,
    _EOL
)

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Join OrthoFinder Orthogroups.tsv files'

num = len


def index_members_by_ortho_ids(ortho):
    index = dict()
    for g in range(num(ortho.groups)):
        for s in range(num(ortho.species)):
            if ortho.groups[g][s] is None:
                continue
            for member in ortho.groups[g][s]:
                index[member] = g
    return index


                


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
    stream.write("  -N,--output-sequence-names\n")
    stream.write("     Map locus names to sequence names internally, then perform filtering.\n")
    stream.write("     Write sequence names to output (use `-n` for locus names).\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized sequences to their\n")
    stream.write("     designated sequence names, not to their placed chromosome names.\n")
    stream.write("     (only effective when --output-sequence-names is enabled)\n")
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
    short_flags = 'hNo:u'
    long_flags = (
        'help',
        'output-sequence-names',
        'output-file=',
        'ignore-unlocalized'
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    is_placed = _placed
    output_seq_names = False
    output_file = sys.stdout
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-o','--output-file'):
            output_file = open(value, 'wt')
        elif flag in ('-N','--output-sequence-names'):
            output_seq_names = True
        elif flag in ('-u','--ignore-unlocalized'):
            is_placed = _localized

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    ortho = OrthoFinderOrthogroups(arguments[0])
    config = SampleConfigFile(arguments[1], load_files=True, map_assigned_molecule=True)
        
    member_indices = index_members_by_ortho_ids(ortho)    
    sample_indices = index_list(ortho.samples)

    singleton_count = 1
    for sample in ortho.samples:
        for locus_name in config.samples[sample.id].loci:
            if locus_name not in member_indices:
                group = ortho.new_group(append=True)
                group[sample.index] = (locus_name,)
                group.id = 'SGL%06d' % singleton_count
                singleton_count += 1

    if output_seq_names:
        for group in ortho.groups:
            for sample in ortho.samples:
                group[sample.index] = sorted(
                    map_loci_to_sequence_counts(
                        group[sample.index],
                        config.samples[sample.id],
                        is_placed,
                        ignore_unplaced=False
                    )
                )

    ortho.to_file(output_file)

    
    
main(sys.argv[1:])
