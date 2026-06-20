
import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Count OrthoFinder Orthogroups.tsv members'


import sys
import getopt
import og.api.count as api

from og.core.members import _LENIENT, _STRICT
from og.core.compression import STDIO, open
from og.core.config import SampleConfigFile
from og.core.assembly_report import is_chr as _localized
from og.core.assembly_report import is_placed as _placed
from og.core.orthogroups.parsers import OrthogroupsFile
from og.constants import _TAB, _EMPTY, _EOL


num = len
_MULTI = '*'  # '\u25CF'


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
    stream.write("     Counts greater than this threshold are converted to the `*` character.\n")
    stream.write("\n")
    stream.write("  -n,--map-to-reference-names\n")
    stream.write("     Map locus names to reference names internally, then calculate counts.\n")
    stream.write("     Requires `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -r,--input-reference-names\n")
    stream.write("     Input orthogroup members are reference names; do not convert, apply\n")
    stream.write("     reference sequence filters appropriately. Requires `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -S,--count-samples\n")
    stream.write("     Instead of member counts, write 1/0 for sample presence/absence\n")
    stream.write("\n")
    stream.write("  -T,--output-totals\n")
    stream.write("     Output counts table with a column of row totals appended.\n")
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
    short_flags = 'hM:o:rSniITuy'
    long_flags = (
        'help',
        'count-samples','count-species',
        'max-count=',
        'output-file=',
        'output-totals',
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

    max_count = sys.maxsize
    config_file = None
    output_file = STDIO
    count_samples = False
    count_totals = False
    ignore_unplaced = 0
    ignore_unlocalized = _placed
    input_reference_names = False
    map_to_reference_names = False
    for flag, value in options:
        if flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-M','--max-count'):
            max_count = int(value)
        elif flag in ('-n','--map-to-reference-names'):
            map_to_reference_names = True
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-r','--input-reference-names'):
            input_reference_names = True
        elif flag in ('-S','--count-samples','--count-species'):
            count_samples = True
        elif flag in ('-T','--output-totals'):
            count_totals = True
        elif flag in ('-u','--ignore-unlocalized'):
            ignore_unlocalized = _localized
        elif flag in ('-y','--yaml'):
            config_file = value

    if map_to_reference_names or input_reference_names:
        if config_file is None:
            usage('`--yaml` config file required')
    else:
        output_reference_names = False
        ignore_unplaced = False            
            
    if num(arguments) != 1:
        usage('Unexpected number of arguments')

    count(
        arguments[0],
        config_file=config_file,
        max_count=max_count,
        count_samples=count_samples,
        count_totals=count_totals,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        map_to_reference_names=map_to_reference_names,
        output_reference_names=output_reference_names,
        output_file=output_file
    )
        
    return 0
    

    
def count(
        ortho_file, config_file=None,
        max_count=sys.maxsize,
        count_samples=False,
        count_totals=False,
        ignore_unplaced=None,
        ignore_unlocalized=False,
        input_reference_names=False,
        map_to_reference_names=False,
        output_reference_names=False,
        output_file=sys.stdout
):
    ortho  = OrthogroupsFile(ortho_file)
    config = None
    if config_file:
        config = SampleConfigFile(
            config_file,
            load_files=True,
            map_assigned_molecule=True
        )
        config.check_samples([s.id for s in ortho.samples])

    counts_table = api.count(
        ortho, config=config,
        count_samples=count_samples,
        input_reference_names=input_reference_names,
        map_to_reference_names=map_to_reference_names
    )

    if output_file:
        output = open(output_file, 'w')
    else:
        output = sys.stdout

    output.write('Orthogroup')
    output.write(_TAB)
    output.write(_TAB.join([s.id for s in counts_table.samples]))
    if count_totals:
        output.write(_TAB)
        output.write('Total')
    output.write(_EOL)
    
    for record in counts_table:
        output.write(record.id)
        for count in record:
            output.write(_TAB)
            output.write(_MULTI if count > max_count else str(count))
        if count_totals:
            output.write(_TAB)
            output.write(str(record.total))
        output.write(_EOL)
            
    output.close()
