#!/usr/bin/env python

import os
import sys

from math import inf
from getopt import getopt, GetoptError
#from fastx.utils.io import load_listfile
from og.core.utils import index_list
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import (
    _COMMA,
    _EMPTY,
    _SPACE,
    _EOL,
    _TAB
)

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Manipulate OrthoFinder orthogroups files'

_LF = '\n'
_CR = '\r'

num = len

def notNone(item):
    return item is not None



def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <qry.tsv> <ref.tsv>\n" % (
        os.path.basename(sys.argv[0])
    ))
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -l,--label-members\n")
    stream.write("     Append the reference orthogroup ID to each query member\n")
    stream.write("\n")
    stream.write("  -o,--output-file <str>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    # stream.write("Notes:\n")
    # stream.write("  - ClusterVenn3 requires `.txt` file suffix\n")
    stream.write("\n\n%s" % message)
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    sys.exit(exitcode)


def main(argv):
    short_flags = 'hlo:'
    long_flags = ('help','label-members','output-file=')
    label_members = False
    output_file = sys.stdout
    try:
        options, arguments = getopt(argv, short_flags, long_flags)
    except GetoptError as message:
        usage(message)

    for flag, value in options:
        if   flag in ('-o','--output-file'):
            output_file = open(value, 'w')
        elif flag in ('-l','--label-members'):
            label_members = True
        elif flag in ('-h','--help'):
            usage(exitcode=0)

    if len(arguments) != 2:
        usage("Unexpected number of arguments")

    qry_ortho = OrthoFinderOrthogroups(arguments[0])
    ref_ortho = OrthoFinderOrthogroups(arguments[1])

    qry_sample_indices = {s.id:s.index for s in qry_ortho.samples}
    ref_sample_indices = {s.id:s.index for s in ref_ortho.samples}

    # 1. Every gene ID in ref_ortho has an Orthogroup ID, create a dict of
    #    these:
    member_to_orthogroup = dict()
    for sample in ref_ortho.samples:
        member_to_orthogroup[sample.id] = dict()
        
        for group in ref_ortho.groups:
            if not group[sample.index]:
                continue

            for member in group[sample.index]:
                member_to_orthogroup[sample.id][member] = group.id
    # 2. Iter through qry_ortho and assign orthogroup IDs given in ref_ortho
    for group in qry_ortho.groups:
        group.id = {group.id: inf}

    for sample in qry_ortho.samples:
        if sample.id not in member_to_orthogroup:
            continue
    
        for group in qry_ortho.groups:
            if not group[sample.index]:
                continue

            members = []
            for member in group[sample.index]:
                ogid = member_to_orthogroup[sample.id].get(member, None)
                
                if ogid is not None:
                    members.append(member + '=' + ogid)
    
                if ogid in group.id:
                    group.id[ogid] += 1
                else:
                    group.id[ogid] = 1

            if label_members:
                group[sample.index] = members

        
    for group in qry_ortho.groups:
        group.id = _COMMA.join(
            filter(notNone, sorted(group.id,key=group.id.get,reverse=True))
        )

    qry_ortho.to_file(output_file)



if __name__ == '__main__':
    main(sys.argv[1:])
