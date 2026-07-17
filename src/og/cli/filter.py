
from og import (
    __authors__,
    __contact__,
    __pkgname__,
    __version__,
    __program__,
)
__purpose__ = 'Filter OrthoFinder Orthogroups.tsv file'

import os
import sys
import getopt
import og.api.filter as api

from og.cli.utils import maxbits, get_exe
from og.core.members import _LENIENT, _STRICT
from og.core.compressio import STDIO
from og.core.config import SampleConfigFile
from og.core.assembly_report import is_chr as _localized
from og.core.assembly_report import is_placed as _placed
from og.core.orthogroups.parsers import OrthogroupsFile
from og.constants import (
    _COMMENT,
    _EMPTY,
    _TAB,
    _EOL
)


num = len

    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    numbits = maxbits()
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
    stream.write("  -m,--min-members <uint>\n")
    stream.write("     Minimum number of members permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -M,--max-members <uint>\n")
    stream.write("     Maximum number of members permitted per orthogroup [2^%d-1]\n" % numbits)
    stream.write("\n")
    stream.write("  -N,--output-reference-names\n")
    stream.write("     Map locus names to reference names internally, then perform filtering.\n")
    stream.write("     Write reference names to output (use `-n` for locus names). Requires\n")
    stream.write("     `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -n,--map-to-reference-names\n")
    stream.write("     Map locus names to reference names internally, then perform filtering.\n")
    stream.write("     Write locus names to output (use `-N` for reference names). Requires\n")
    stream.write("     `--yaml` be defined.\n")
    stream.write("\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -r,--input-reference-names\n")
    stream.write("     Input orthogroup members are reference names; do not convert, apply\n")
    stream.write("     reference sequence filters appropriately. Requires `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -s,--min-samples <uint>\n")
    stream.write("     Minimum number of samples permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -S,--max-samples <uint>\n")
    stream.write("     Maximum number of samples permitted per orthogroup [2^%d-1]\n" % numbits)
    stream.write("\n")
    stream.write("  -t,--tree-filter\n")
    stream.write("     Filter orthogroups using the in.yaml ploidy tree.\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized references to their\n")
    stream.write("     original reference names, not to their assigned molecule sequence names.\n")
    stream.write("\n")
    stream.write("  -v,--invert-output\n")
    stream.write("     Output orthogroups failing the specified filters [passing]\n")
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
    stream.write("  2. The ploidy submapping of the `tree` key in in.yaml is a Newick-formatted\n")
    stream.write("     tree string that can be written to filter orthogroups by conditioning\n")
    stream.write("     on numbers of members (loci or references) and samples at each leaf\n")
    stream.write("     node and internal node, respectively. In place of Newick branch lengths,\n")
    stream.write("     however, the admissible number of members and samples are specified\n")
    stream.write("     using unsigned integer number ranges, consisting (inclusively) of the\n")
    stream.write("     minimum number, a dash (`-`), then maximum number without any intervening\n")
    stream.write("     whitespace. The minimum or maximum may be omitted to specify open ranges.\n")
    stream.write("     For example, the tree below can be interpreted as follows:\n")
    stream.write("         '((A:1, B:-4):1-2, (C:0-1, D):1-):2-4'\n")
    stream.write("     Leaf node A must have one, and exactly one, member present; leaf node\n")
    stream.write("     B may have up to four members (inclusive); the A+B subclade (here,\n")
    stream.write("     represented as an internal node) requires one-to-two (inclusive)\n")
    stream.write("     samples to be be present. Leaf node C must be present at most once.\n")
    stream.write("     The number of members on leaf node D is unrestricted; and subclade C+D\n")
    stream.write("     requires at least one samples be present. Because of the two subclade-\n")
    stream.write("     specific constraints, the two-to-four samples required at the root\n")
    stream.write("     will always also be satisfied.\n")
    stream.write("\n")    
    # #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    # #            0        10        20        30        40        50        60        70        80
    # stream.write("  3. A locus BED table is a two-column file specifying the samples names\n")
    # stream.write("     (same as used in the Orthogroups.tsv header) and paths to locus BED\n")
    # stream.write("     files. The BED files must contain locus names in fourth column and the\n")
    # stream.write("     samples names prepended to the reference names (e.g., Hsa1 for chromosomes\n")
    # stream.write("     and HsaSca123 or HsaUn123 for unplaced scaffolds). If a locus BED table\n")
    # stream.write("     is given, then the number of references per samples is counted as the\n")
    # stream.write("     members rather than locus names.\n")
    # stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    

def main(argv):
    short_flags = 'hIio:m:M:Nnrs:S:tuvy:'
    long_flags = (
        'help',
        'output-file=',
        'ignore-unplaced-leniently',
        'ignore-unplaced-strictly',
        'ignore-unlocalized',
        'ignore-ploidy-tree',
        'input-reference-names',
        'output-reference-names',
        'map-to-reference-names',
        'min-members=',
        'max-members=',
        'min-samples=','min-species=',
        'max-samples=','max-species=',
        'invert-output',
        'yaml='
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    invert = False
    config_file = None
    output_file = STDIO
    min_members = 1
    max_members = sys.maxsize
    min_samples = 1
    max_samples = sys.maxsize
    map_to_reference_names = False
    input_reference_names = False
    output_reference_names = False
    ignore_unplaced = 0
    ignore_unlocalized = _placed
    ploidy_tree_filter = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-m','--min-members'):
            min_members = int(value)
        elif flag in ('-M','--max-members'):
            max_members = int(value)
        elif flag in ('-N','--output-reference-names'):
            output_reference_names = True
            map_to_reference_names = True
        elif flag in ('-n','--map-to-reference-names'):
            map_to_reference_names = True
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-s','--min-samples','--min-species'):
            min_samples = int(value)
        elif flag in ('-S','--max-samples','--max-species'):
            max_samples = int(value)
        elif flag in ('-r','--input-reference-names'):
            input_reference_names = True
        elif flag in ('-t','--tree-filter'):
            ploidy_tree_filter = True
        elif flag in ('-u','--ignore-unlocalized'):
            ignore_unlocalized = _localized            
        elif flag in ('-v','--invert-output'):
            invert = True
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


    filter(
        arguments[0],
        config_file=config_file,
        invert=invert,
        min_samples=min_samples,
        max_samples=max_samples,
        min_members=min_members,
        max_members=max_members,
        ploidy_tree_filter=ploidy_tree_filter,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        map_to_reference_names=map_to_reference_names,
        output_reference_names=output_reference_names,
        output_file=output_file
    )

    return 0
    
    

def filter(
        ortho_file,
        config_file=None,
        invert=False,
        min_samples=1,
        max_samples=sys.maxsize,
        min_members=1,
        max_members=sys.maxsize,
        ploidy_tree_filter=False,
        ignore_unplaced=False,
        ignore_unlocalized=False,
        input_reference_names=False,
        map_to_reference_names=False,
        output_reference_names=False,
        output_file=sys.stdout
):
    ortho = OrthogroupsFile(ortho_file)    
    config = None
    if config_file is None:
        ploidy_tree_filter = False
    else:
        config = SampleConfigFile(
            config_file,
            load_files=True,
            map_assigned_molecule=True
        )
        config.check_samples([s.id for s in ortho.samples])

    output = api.filter(
        ortho, config=config,
        invert=invert,
        min_samples=min_samples,
        max_samples=max_samples,
        min_members=min_members,
        max_members=max_members,
        ploidy_tree_filter=ploidy_tree_filter,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        map_to_reference_names=map_to_reference_names,
        output_reference_names=output_reference_names
    )

    output.to_file(output_file)
