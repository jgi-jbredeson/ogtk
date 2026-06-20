#!/usr/bin/env python


import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Reassign groups using posterior probabilities'


import io
import sys
import getopt
import og.api.update as api

from og.cli.utils import read_list
from og.core.compression import STDIO, open
from og.core.orthogroups.parsers import CountedClusteredOrthogroups
from og.core.orthogroups.parsers import ClusterErrorOrthogroups
from og.constants import (
    _COMMA,
    _COMMENT,
    _EMPTY,
    _EOL,
    _SPACE,
    _TAB,
)


num = len


def format_grouptag(ortho, group):
    return '%s%s' % (
        ortho.grouptag,
        -1 if group.cluster is None else group.cluster
    )
    
    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <training.tsv> <error.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -G,--ignore-group-names <name1>[,<name2>[,...]]\n")
    stream.write("     Exclude groups, designated by their comma-separated list of group\n")
    stream.write("     names, from reassignment. Orthogroups in the designated synteny groups\n")
    stream.write("     will not be reassigned.\n")
    stream.write("\n")
    stream.write("  -p,--min-prob-reassign <float>  (default: 0.5)\n")
    stream.write("     Reassign orthogroups with posterior probabilities greater-than or equal\n")
    stream.write("     to the specified threshold.\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  training.tsv is the clustered input file passed to score-orthogroups, with\n")
    stream.write("  lines belonging to the same group sorted together and the group preceded\n")
    stream.write("  by `##group=N` meta lines to demarcate the groups and define their cluster\n")
    stream.write("  IDs (N).\n")
    stream.write("\n")
    stream.write("  error.tsv is the output of score-orthogroups run with the -E/--check-errors\n")
    stream.write("  flag enabled.\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)



def main(argv):
    short_options = 'hG:p:'
    long_options = (
        'help',
        'ignore-group-names=',
        'min-prob-reassign='
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    output_file = STDIO
    ignore_groups = set()
    min_prob_reassign = 0.5
    min_prob_preserve = -1.0
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-o','--output-file'):
            output_file = value
        elif flag in ('-G','--ignore-group-names'):
            ignore_groups = set(read_list(value))
        elif flag in ('-p','--min-prob-reassign'):
            min_prob_reassign = float(value)

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    update(
        training_file,
        error_file,
        min_prob_reassign=min_prob_reassign,
        ignore_groups=ignore_groups
    )

    return 0



def update(training_file, error_file):
    training_ortho = OrthogroupsFile(training_file)
    error_ortho = ClusterErrorOrthogroups(error_file)
    
    updated_ortho = api.update(training_ortho, error_ortho)

    updated_ortho.groups.sort(key=lambda g: (g.cluster, -g.count))
    with open(output_file, 'w') as output:
        print(updated_ortho.format_header())
        cluster = None
        for group in updated_ortho:
            if group.cluster != cluster:
                print(format_grouptag(updated_ortho, group), file=output)
            cluster = group.cluster
            print(updated_ortho.format_record(group=group), file=output)
