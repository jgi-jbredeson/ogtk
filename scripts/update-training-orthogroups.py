#!/usr/bin/env python

# TODO: modify to allow reading different file input Orthofinder-like file formats

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

from og.core.parsers.orthogroups import CountedClusteredOrthogroups
from og.core.parsers.orthogroups import ClusterErrorOrthogroups
from og.constants import (
    _COMMA,
    _COMMENT,
    _EMPTY,
    _EOL,
    _SPACE,
    _TAB,
)

_DEBUG = False

num = len



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
    stream.write("  -G,--ignore-groups <lab1>[,<lab2>[,...]]\n")
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
        'ignore-groups=',
        'min-prob-reassign='
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    ignore_groups = set()
    min_prob_reassign = 0.5
    min_prob_preserve = -1.0
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-G','--ignore-groups'):
            ignore_groups = set(value.rstrip(_COMMA).split(_COMMA))
        elif flag in ('-p','--min-prob-reassign'):
            min_prob_reassign = float(value)
        # elif flag in ('-P','--max-prob-preserve'):
        #     min_prob_preserve = float(value)

    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    orthoM = CountedClusteredOrthogroups(arguments[0])
    orthoP = ClusterErrorOrthogroups(arguments[1])

    discard = dict()
    reassign = dict()
    preserve = dict()
    for group in orthoP.groups:
        if group.cluster[1] in ignore_groups:
            continue
        if group.probability[1] >= min_prob_reassign:
            print(group.id)
            reassign[group.id] = (group.cluster[1], group.probability[1])

    sortorder = dict()
    for g, group in enumerate(orthoM.groups):
        if group.id in reassign:
            group.cluster = reassign[group.id][0]
            group.probability = reassign[group.id][1]

        if group.cluster not in sortorder:
            sortorder[group.cluster] = []
            
        sortorder[group.cluster].append((-1*int(group.count), g))

    print(orthoM.format_orthogroups_header())
    for group_id in sortorder:
        sortorder[group_id].sort()

        print('%s%s' % (orthoM.grouptag, group_id))
        for count, g in sortorder[group_id]:
            print(orthoM.format_orthogroups_record(index=g))
            

main(sys.argv[1:])
