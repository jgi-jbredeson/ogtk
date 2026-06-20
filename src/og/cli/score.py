
import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Use manual clustering to classify new groups'


import sys
import getopt
import og.api.score as api

from og.constants import _COMMA, _EMPTY
from og.cli.utils import read_list
from og.core.config import SampleConfigFile
from og.core.members import _LENIENT, _STRICT
from og.core.compression import STDIO
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
    stream.write("Usage: %s [options] <classified.tsv> [unclassified.tsv] <in.yaml>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80    
    stream.write("  -E,--check-errors\n")
    stream.write("     Check for classification errors in classified.tsv\n")
    stream.write("\n")
    stream.write("  -g,--ignore-group-names <name1>[,<name2>[,...]]\n")
    stream.write("     Exclude groups, designated by their comma-separated list of group\n")
    stream.write("     names, from forming Bayesian classification groups. Orthogroups in\n")
    stream.write("     the designated synteny groups will be clustered to other groups.\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     When mapping locus names to reference names, remove unplaced reference\n")
    stream.write("     members. See also `--ignore-unlocalized`.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     When mapping locus names to reference names, convert unplaced reference\n")
    stream.write("     names to the value of `unplaced_id` in the YAML config for that species.\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -P,--output-posterior-prob-file <file>\n")
    stream.write("     Write the Bayesian posterior probabilities table to file.\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized reference to their\n")
    stream.write("     original reference names, not to their assigned molecule sequence names.\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  classified.tsv is a mrg.counts.tsv file output by cluster-orthogroups.py\n")
    stream.write("  or an OrthoFinder orthogroups.tsv-formatted file. Lines belonging to the\n")
    stream.write("  same syntenic groups are assumed to be grouped together and such groups\n")
    stream.write("  preceded by `##group=N` meta lines to demarcate the groups and define their\n")
    stream.write("  group/cluster names (in this example, `N`).\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    

def main(argv):
    short_options = 'hEg:Iio:P:u'
    long_options = (
        'help',
        'check-errors',
        'ignore-group-names=',
        'ignore-unlocalized',
        'ignore-unplaced-strictly',
        'ignore-unplaced-liniently',
        'output-posterior-prob-file=',
        'output-file=',
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    check_errors = False
    output_file = STDIO
    ignore_groups = set()
    ignore_unplaced = False
    ignore_unlocalized = _placed
    output_posterior_file = None
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-E','--check-errors'):
            check_errors = True
        elif flag in ('-g','--ignore-group-names'):
            ignore_groups = set(read_list(value))
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-i','--ignore-unplaced-liniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-P','--output-posterior-prob-file'):
            output_posterior_file = value
        elif flag in ('-u','--ignore-unlocalized'):
            ignore_unlocalized = _localized
            
    if num(arguments) != 2 and \
       num(arguments) != 3:
        usage('Unexpected number of arguments')

    if num(arguments) == 2:
        arguments.append(arguments[1])
        arguments[1] = arguments[0]
    if check_errors:
        arguments[1] = arguments[0]

    score(
        classify_file=arguments[1],
        training_file=arguments[0],
        config_file=arguments[2],
        ignore_groups=ignore_groups,
        check_errors=check_errors,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        output_posterior_file=output_posterior_file,
        output_file=output_file
    )

    return 0



def score(
        classify_file, training_file, config_file,
        check_errors=False, ignore_groups=(),
        ignore_unplaced=False, ignore_unlocalized=_placed,
        output_posterior_file=None, output_file=sys.stdout
):
    if classify_file is None:
        classify_file = training_file
    if check_errors:
        classify_file = training_file

    training_ortho = OrthogroupsFile(training_file)
    classify_ortho = OrthogroupsFile(classify_file)
    config = SampleConfigFile(config_file,
                              load_files=True, map_assigned_molecule=True)
    
    scored, pprobs = \
        api.score(
            classify_ortho,
            training_ortho,
            config,
            check_errors=check_errors,
            ignore_groups=ignore_groups,
            ignore_unplaced=ignore_unplaced,
            ignore_unlocalized=ignore_unlocalized
        )
    
    if output_posterior_file:
        pprobs.to_file(output_posterior_file)
    scored.to_file(output_file)

