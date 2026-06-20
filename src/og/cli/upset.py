
import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Print upset plot of Orthogroups.tsv file'


import sys
import getopt
import og.api.upset as api

from og.core.members import _LENIENT, _STRICT
from og.core.compression import STDIO, open
from og.core.config import SampleConfigFile

from og.core.assembly_report import is_chr as _localized
from og.core.assembly_report import is_placed as _placed
from og.core.orthogroups.parsers import OrthogroupsFile
from og.core.orthogroups.formatters import (
    GROUP_MEMBERSHIP,
    GROUP_MULTIPLES,
    SORT_CLUSTERS,
    SORT_MEMBERS
)
from og.constants import (
    _EMPTY,
    _SPACE,
    _EOL,
    _TAB
)


num = len


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <in.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -a,--force-ascii\n")
    stream.write("     Force writing output in ASCII-only characters\n")
    stream.write("\n")
    stream.write("  -g,--group-by <enum>\n")
    stream.write("     Group output by (1) membership, (2) number of multiples\n")
    stream.write("     (See the `--max-count` option below), or perform no grouping [0]\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     Strictly ignore unplaced references in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains no chromosomal references, that cell\n")
    stream.write("     then contains no members.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     Leniently ignore unplaced references in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains only unplaced (ie, non-chomosomal)\n")
    stream.write("     references, that cell contains members.\n")
    stream.write("\n")
    stream.write("  -M,--max-count <uint>\n")
    stream.write("     Count up to `-M` number of references per sample and orthogroup.\n")
    stream.write("     Counts greater than this threshold are converted to the `\u25CF` character\n")
    stream.write("     (or `*` if `--force-ascii` is enabled).\n")
    stream.write("\n")
    stream.write("  -n,--map-to-reference-names\n")
    stream.write("     Map locus names to reference names internally, then calculate plot.\n")
    stream.write("     Requires `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -r,--input-reference-names\n")
    stream.write("     Input orthogroup members are reference names; do not convert, apply\n")
    stream.write("     reference sequence filters appropriately. Requires `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -s,--sort-by <enum>\n")
    stream.write("     Sort plot rows by (1) number of clusters or (2) number of members [-1]\n")
    stream.write("     (Negate the enumerative value to reverse sort order)\n")
    stream.write("\n")
    stream.write("  -S,--upset-separator <str>\n")
    stream.write("     Add space between columns of the upset plot using the given separator\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized references to their\n")
    stream.write("     original reference names, not to their assigned molecule sequence names.\n")
    stream.write("\n")
    stream.write("  -y,--yaml <in.yaml>\n")
    stream.write("     Input YAML config file, with `tree` mapping key with nested submapping\n")
    stream.write("     `ploidy` key to a Newick tree string value. For each species, a mapping\n")
    stream.write("     with required submappings: `loci` and `references` keys with file path\n")
    stream.write("     values, and `unplaced_id` key with string value.\n")
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
    short_flags = 'ahS:s:g:M:o:rniIuy:'
    long_flags = (
        'force-ascii',
        'help',
        'upset-sep=',
        'sort-by=',
        'group-by=',
        'max-count=',
        'output-file=',
        'input-reference-names',
        'map-to-reference-names',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'ignore-unlocalized',
        'yaml='
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    sort_by = -1
    group_by = 0
    max_count = 0
    force_ascii = False
    output_file = STDIO
    config_file = None
    upset_separator = _EMPTY
    ignore_unplaced = 0
    ignore_unlocalized = _placed
    map_to_reference_names = False
    input_reference_names = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-a','--force-ascii'):
            upset_separator = _SPACE
            force_ascii = True
        elif flag in ('-g','--group-by'):
            group_by |= (0x1 | int(value))
        elif flag in ('-s','--sort-by'):
            sort_by = int(value)
        elif flag in ('-S','--upset-separator'):
            upset_separator = value.encode('utf-8').decode('unicode_escape')
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-M','--max-count'):
            max_count = int(value)
        elif flag in ('-n','--map-to-reference-names'):
            map_to_reference_names = True
        elif flag in ('-r','--input-reference-names'):
            input_reference_names = True
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-u','--ignore-unlocalized'):
            ignore_unlocalized = _localized
        elif flag in ('-y','--yaml'):
            config_file = value

    if not (abs(sort_by) & (SORT_CLUSTERS | SORT_MEMBERS)):
        usage("Invalid enumerative value: --sort-by=%s" % str(sort_by))
        
    if map_to_reference_names or input_reference_names:
        if config_file is None:
            usage('`--yaml` config file required')
    else:
        ignore_unplaced = False            
            
    if num(arguments) != 1:
        usage('Unexpected number of arguments')

    upset(
        arguments[0],
        config_file=config_file,
        max_count=max_count,
        sort_by=sort_by,
        group_by=group_by,
        force_ascii=force_ascii,
        upset_separator=upset_separator,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        map_to_reference_names=map_to_reference_names,
        output_file=output_file
    )

    return 0
    


def upset(
        ortho_file,
        config_file=None,
        max_count=0,
        sort_by=0,
        group_by=0,
        force_ascii=False,
        upset_separator=_EMPTY,
        ignore_unplaced=False,
        ignore_unlocalized=None,
        input_reference_names=False,
        map_to_reference_names=False,
        output_file=sys.stdout
):
    ortho = OrthogroupsFile(ortho_file)
    config = None
    if config_file:
        config = SampleConfigFile(
            config_file,
            load_files=True,
            map_assigned_molecule=True
        )

    upset = api.upset(
        ortho,
        config=config,
        max_count=max_count,
        sort_by=sort_by,
        group_by=group_by,
        force_ascii=force_ascii,
        upset_separator=upset_separator,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        map_to_reference_names=map_to_reference_names
    )

    with open(output_file, 'w') as output:
        output.write(upset + _EOL)

