#!/usr/bin/env python

import os
import sys

from getopt import getopt, GetoptError
from og.core.utils import index_list
from og.core.compression import is_stream
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
_SORT_OGIDS = 0x1
_SORT_MEMBERS = 0x2
_SORT_SAMPLES = 0x4

_TWO_SAMPLES_REQUIRED = \
    "Exactly two samples required with `--output-type A`"
_OGID_FIELD = {'Orthogroup','OG','HOG'}


num = len


def lineparser(line):
    return line.rstrip(_EOL).split(_TAB)[0]


def load_listfile(listfiles, parser=lambda x: x):
    """Load a non-redundant list of items from listfiles.

    When multiple instances of the same value item are encountered,
    just the first is kept.

    parser is a callable that inputs a single string and returns 
    an object to be appended to the end of the growing list.
    """
    if isinstance(listfiles, str):
        listfiles = (listfiles,)
    itemset = set()
    itemlist = list()
    for listfile in listfiles:
        inputfile = open(listfile, 'r')
        for line in inputfile:
            item = parser(line)

            if item in itemset:
                continue
        
            itemlist.append(item)
            itemset.add(item)

        inputfile.close()

    return itemlist


def read_order_file(filename):
    return load_listfile(filename, parser=lineparser)
    # if "Orthogroup" not in order:
    #     order = ["Orthogroup"] + order
    # return order


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <Orthogroups.tsv>\n" % (
        os.path.basename(sys.argv[0])
    ))
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -d,--prefix-delim <char>\n")
    stream.write("     When prefixing sample names to locus IDs, seperate them using char [|]\n")
    stream.write("\n")
    stream.write("  -D,--replace-delim <char>\n")
    stream.write("     To make IDs safe, replace existing specified char in locus IDs prior to\n")
    stream.write("     prefixing sample ID [|]\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -O,--output-type <str>\n")
    stream.write("     Output orthogroups in specified format: [F], OrthoFinder; V, OrthoVenn;\n")
    stream.write("     A, MCScan anchors format\n")
    stream.write("\n")
    stream.write("  -p,--prefix-sample-names\n")
    stream.write("     Prepend the locus IDs with the sample names declared in the header.\n")
    stream.write("\n")
    stream.write("  -P,--remove-sample-names\n")
    stream.write("     Remove prefixed sample name {-p} and delimiter {-d} from locus IDs\n")
    stream.write("\n")
    stream.write("  -s,--sample-order <sample1[,sample2[,...]]>\n")
    stream.write("     Input comma-separated list of desired output sample order.\n")
    stream.write("\n")    
    stream.write("  -S,--sample-order-file <file>\n")
    stream.write("     Input file listing (one per line) the desired output sample order.\n")
    stream.write("\n")
    stream.write("  --sort-by <enum>\n")
    stream.write("     Sort rows by (1) orthogroup IDs, numbers of (2) members or (3) samples,\n")
    stream.write("     or maintain input sort order [0].\n")
    stream.write("     (Negate the enumerative value to reverse sort order)\n")    
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this usage message\n")
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  - ClusterVenn3 requires `.txt` file suffix\n")
    stream.write("\n\n%s" % message)
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    sys.exit(exitcode)
    
    
def main(argv):
    err = sys.stderr
    sort_by = 0
    output_file = sys.stdout
    prefix_delim = '|'
    replace_delim = prefix_delim
    prefix_samples = False
    remove_samples = False
    oldsamples = []
    oldindices = {}
    newsamples = []
    newindices = {}
    output_type = 'F'
    valid_output_types = set('AFV')
    long_flags = (
        'outfile=',
        'output-file=',
        'output-type=',
        'sample-order=','species-order=',
        'sample-order-file=','species-order-file=',
        'replace-delim',
        'prefix-delim=',
        'prefix-sample-names','prefix-species-names',
        'remove-sample-names','remove-species-names',
        'sort-by=',
        'orthovenn',
        'help'
    )
    short_flags = 'd:o:O:D:PpS:s:vh'
    
    try:
        options, arguments = getopt(argv, short_flags, long_flags)
    except GetoptError as message:
        usage(message)

    for flag, value in options:
        if   flag in ('-o','--outfile','--output-file'):
            output_file = open(value, 'w')
        elif flag in ('-O','--output-type'):
            output_type = value
        elif flag in ('-d','--prefix-delim'):
            prefix_delim = value
        elif flag in ('-D','--replace-delim'):
            replace_delim = value
        elif flag in ('-p','--prefix-sample-names','--prefix-species-names'):
            prefix_samples = True
            remove_samples = False
        elif flag in ('-P','--remove-sample-names','--remove-species-names'):
            prefix_samples = False
            remove_samples = True
        elif flag in ('-S','--sample-order-file','--species-order-file'):
            newsamples.extend(read_order_file(value))
        elif flag in ('-s','--sample-order','--species-order'):
            newsamples.extend(value.split(_COMMA))
        elif flag in ('--sort-by',):
            sort_by = int(value)
        elif flag in ('-v','--orthovenn'):
            output_type = 'V'
        elif flag in ('-h','--help'):
            usage(exitcode=0)

    if output_type not in valid_output_types:
        usage("Unsupported output type: `%s`" % output_type)
    if sort_by:
        if not (abs(sort_by) & (_SORT_OGIDS | _SORT_MEMBERS | _SORT_SAMPLES)):
            usage("Invalid enumerative value: --sort-by=%s" % str(sort_by))
        
    if len(arguments) != 1:
        usage("Unexpected number of arguments")
        
    ortho = OrthoFinderOrthogroups(arguments[0])
    
    oldsamples = [s.id for s in ortho.samples]
    oldindices = {s.id:s.index for s in ortho.samples}

    ogid_as_sample = False
    if newsamples:
        _newsamples = []
        for sample in newsamples:
            if sample in _OGID_FIELD:
                ogid_as_sample = True
                continue
            if sample not in oldindices:
                raise KeyError(
                    "Sample not found in orthogroups "
                    "file: `%s`" % sample
                )
            _newsamples.append(sample)
        
        newsamples = _newsamples
        newindices = index_list(newsamples)
    else:
        newsamples = oldsamples
        newindices = oldindices

    _ortho = OrthoFinderOrthogroups(samples=newsamples)
    for oldgroup in ortho.groups:
        newgroup = _ortho.new_group(append=False)
        newgroup.id = oldgroup.id
        passes = False
        for sample in newsamples:
            if num(oldgroup[oldindices[sample]]) > 0:
                newgroup[newindices[sample]] = oldgroup[oldindices[sample]]
                passes = True
        if passes:
            _ortho.groups.append(newgroup)
    ortho = _ortho

    if prefix_samples:
        for group in ortho.groups:
            for sample in ortho.samples:
                if not group[sample.index]:
                    continue
                members = list()
                prefix = sample.id + prefix_delim
                for member in group[sample.index]:
                    if not member.startswith(prefix):
                        member = prefix + member.replace(prefix_delim,replace_delim)
                    members.append(member)
                group[sample.index] = members
                
    if remove_samples:
        for group in ortho.groups:
            for sample in ortho.samples:
                if not group[sample.index]:
                    continue
                members = list()
                prefix = sample.id + prefix_delim
                for member in group[sample.index]:
                    if member.startswith(prefix):
                        member = member[len(prefix):]
                    members.append(member)
                group[sample.index] = members
            
    if output_type == 'V':
        ortho.groups.sort(key=lambda g: -g.num_members)
        for group in ortho.groups:
            sep = _EMPTY
            if ogid_as_sample:
                output_file.write(group.id)
                sep = _TAB
            for sample in ortho.samples:
                if group[sample.index]:
                    output_file.write(
                        sep + _TAB.join(group[sample.index])
                    )
                    sep = _TAB
            output_file.write(_LF)

    elif output_type == 'A':
        if ogid_as_sample:
            if num(ortho.samples) != 1:
                raise Exception(_TWO_SAMPLES_REQUIRED)
        elif num(ortho.samples) != 2:
            raise Exception(_TWO_SAMPLES_REQUIRED)

        for group in ortho.groups:
            if ogid_as_sample:
                sample_i = (group.id,)
                sample_j = group[0]
            else:
                sample_i = group[0]
                sample_j = group[1]

            for member_i in sample_i:
                for member_j in sample_j:
                    output_file.write(
                        _TAB.join((member_i, member_j)) + _LF
                    )
    else:
        if sort_by:
            if abs(sort_by) & _SORT_OGIDS:
                key = lambda g: g.id
            elif abs(sort_by) & _SORT_MEMBERS:
                key = lambda g: g.num_members
            elif abs(sort_by) & _SORT_SAMPLES:
                key = lambda g: g.num_samples
            else:
                raise AssertionError('Invalid sort value: %s' % str(sort_by))
            ortho.groups.sort(key=key, reverse=(sort_by < 0))
                              
        ortho.to_file(output_file)

    if is_stream(output_file):
        output_file.close()
        
main(sys.argv[1:])
