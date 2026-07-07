
from og import (
    __authors__,
    __contact__,
    __pkgname__,
    __version__,
    __program__,
)
__purpose__ = 'Append missing members to an orthogroups file'


import os
import sys
import getopt
import og.api.add_singletons as api

from og.constants import _EMPTY
from og.cli.utils import get_exe
from og.core.compressio import STDIO
from og.core.config import SampleConfigFile
from og.core.assembly_report import is_chr as _localized
from og.core.assembly_report import is_placed as _placed
from og.core.orthogroups.parsers import OrthogroupsFile



num = len



def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] --yaml <in.yaml> <in.tsv>\n" % get_exe(__program__, __file__))
    stream.write("\n")
    stream.write("Required:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -y,--yaml <in.yaml>\n")
    stream.write("     Input YAML config file, with `tree` mapping key with nested submapping\n")
    stream.write("     `ploidy` key to a Newick tree string value. For each species, a mapping\n")
    stream.write("     with required submappings: `loci` and `references` keys with file path\n")
    stream.write("     values, and `unplaced_id` key with string value.\n")
    stream.write("\n")    
    stream.write("Options:\n")
    stream.write("  -N,--output-reference-names\n")
    stream.write("     Map locus names to reference names internally, then perform filtering.\n")
    stream.write("     Write reference names to output (use `-n` for locus names).\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized references to their\n")
    stream.write("     original reference names, not to their assigned molecule sequence names. (only\n")
    stream.write("     effective when --output-reference-names is enabled)\n")
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
    short_flags = 'hNo:uy:'
    long_flags = (
        'help',
        'output-reference-names',
        'output-file=',
        'ignore-unlocalized',
        'yaml='
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    config_file = None        
    output_file = STDIO
    ignore_unlocalized = _placed
    output_reference_names = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-N','--output-reference-names'):
            output_reference_names = True
        elif flag in ('-u','--ignore-unlocalized'):
            ignore_unlocalized = _localized
        elif flag in ('-y','--yaml'):
            config_file = value

    if config_file is None:
        usage('`--yaml` config file required')
    if num(arguments) != 1:
        usage('Unexpected number of arguments')

    add_singletons(
        arguments[0],
        config_file,
        ignore_unlocalized=ignore_unlocalized,
        output_reference_names=output_reference_names,
        output_file=output_file
    )

    return 0


def add_singletons(
        ortho_file, config_file,
        ignore_unlocalized=None,
        output_reference_names=False,
        output_file=sys.stdout
):
    ortho = OrthogroupsFile(ortho_file)
    config = SampleConfigFile(
        config_file,
        load_files=True,
        map_assigned_molecule=True
    )

    api.add_singletons(
        ortho, config,
        output_reference_names=output_reference_names,
        ignore_unlocalized=ignore_unlocalized
    )
    
    ortho.to_file(output_file)

