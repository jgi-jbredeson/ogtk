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
    short_options = 'hP:'
    long_options = (
        'help',
        'min-prob-reassign='
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    min_prob_reassign = 0.5
    min_prob_preserve = -1.0
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-P','--min-prob-reassign'):
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
    for i in range(num(orthoP.ids)):
        if orthoP.probabilities[i][1] >= min_prob_reassign:
            reassign[orthoP.ids[i]] = (orthoP.clusters[i][1], orthoP.probabilities[i][1])
        # elif orthoP.probabilities[i][1] >= min_prob_preserve:
        #     preserve[orthoP.ids[i]] = orthoP.clusters[0]

    sortorder = dict()
    for i in range(num(orthoM.ids)):
        if orthoM.ids[i] in reassign:
            orthoM.clusters[i] = reassign[orthoM.ids[i]][0]
            orthoM.probabilities[i] = reassign[orthoM.ids[i]][1]

        if orthoM.clusters[i] not in sortorder:
            sortorder[orthoM.clusters[i]] = []
            
        sortorder[orthoM.clusters[i]].append((-1*int(orthoM.counts[i]), i))

    print(orthoM.format_orthogroups_header())
    for group_id in sortorder:
        sortorder[group_id].sort()

        print('%s%s' % (orthoM.grouptag, group_id))
        for count, i in sortorder[group_id]:
            print(orthoM.format_orthogroups_record(index=i))
            

main(sys.argv[1:])
