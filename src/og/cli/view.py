#!/usr/bin/env python

import os
import sys
import og.api.view as api

from getopt import getopt, GetoptError
from og.core.compression import STDIO, open
from og.core.members import _LENIENT, _STRICT
from og.core.tsv import read_to_list
from og.core.config import SampleConfigFile
from og.core.assembly_report import is_chr as _localized
from og.core.assembly_report import is_placed as _placed
from og.core.orthogroups.parsers import OrthogroupsFile
from og.core.orthogroups.formatters import OrthoVennFormatter
from og.core.orthogroups.formatters import MCScanAnchorsFormatter
from og.api.view import ENUM_OUTPUT, ENUM_SORTBY, ENUM_CALLKEY
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

_OGID_FIELD = {'Orthogroup','OG','HOG'}
_TWO_SAMPLES_REQUIRED = \
    "Exactly two samples required with `--output-type A`"


num = len

def read_listfiles(listfiles):
    """Load a non-redundant list of items from listfiles.

    When multiple instances of the same value item are encountered,
    just the first is kept.

    parser is a callable that inputs a single string and returns 
    an object to be appended to the end of the growing list.
    """
    if isinstance(listfiles, (bytes, str)):
        listfiles = (listfiles,)
    itemset = set()
    itemlist = list()
    for listfile in listfiles:
        for item in read_to_list(listfile):
            if item in itemset:
                continue
            itemlist.append(item)
            itemset.add(item)

    return itemlist



def has_ogid(sample_names):
    for sample_name in sample_names:
        if sample_name in _OGID_FIELD:
            return sample_name
    return False



def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s [options] <in.tsv>\n" % (
        os.path.basename(sys.argv[0])
    ))
    stream.write("\n")
    stream.write("Options:\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80    
    stream.write("  -D,--replace-delim <char>\n")
    stream.write("     To make IDs safe, replace existing specified char in locus IDs prior to\n")
    stream.write("     prefixing sample ID [|]\n")
    stream.write("\n")
    stream.write("  -d,--prefix-delim <char>\n")
    stream.write("     When prefixing sample names to locus IDs, seperate them using char [|]\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     When mapping locus names to reference names, remove unplaced reference\n")
    stream.write("     members. See also `--ignore-unlocalized`.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     When mapping locus names to reference names, convert unplaced reference\n")
    stream.write("     names to the value of `unplaced_id` in the YAML config for that species.\n")
    stream.write("\n")
    stream.write("  -N,--output-reference-names\n")
    stream.write("     Map locus names to reference names and write to output. Requires `--yaml`\n")
    stream.write("     be defined.\n")
    stream.write("\n")
    stream.write("  -O,--output-type <str>\n")
    stream.write("     Output orthogroups in specified format: [F], OrthoFinder; V, OrthoVenn;\n")
    stream.write("     A, MCScan anchors format\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -P,--remove-sample-prefix\n")
    stream.write("     Remove prefixed sample names {-p} and delimiter {-d} from member IDs\n")
    stream.write("\n")
    stream.write("  -p,--prepend-sample-prefix\n")
    stream.write("     Prepend sample name and delimiter {-d} to each member ID. Sample names\n")
    stream.write("     are determined from the in.tsv file header.\n")
    stream.write("\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("  -r,--input-reference-names\n")
    stream.write("     Input orthogroup members are reference names; do not convert, apply\n") 
    stream.write("     reference sequence filters appropriately. Requires `--yaml` be defined.\n")
    stream.write("\n")
    stream.write("  -S,--samples-file <file>\n")
    stream.write("     Input file listing (one per line) the samples to output (ordered).\n")
    stream.write("\n")
    stream.write("  -s,--samples-list <sample1[,sample2[,...]]>\n")
    stream.write("     Input comma-separated list of samples to output (ordered).\n")
    stream.write("\n")    
    stream.write("  --sort-by <enum>\n")
    stream.write("     Sort rows by (1) orthogroup IDs, numbers of (2) members or (4) samples,\n")
    stream.write("     or maintain input sort order [0].\n")
    stream.write("     (Negate the enumerative value to reverse sort order)\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized reference to their\n")
    stream.write("     original reference names, not to their assigned molecule sequence names.\n")
    stream.write("\n")
    stream.write("  -y,--yaml <in.yaml>\n")
    stream.write("     Input YAML config file, with `tree` mapping key with nested submapping\n")
    stream.write("     `ploidy` key to a Newick tree string value. For each species, a mapping\n")
    stream.write("     with required submappings: `loci` and `references` keys with file path\n")
    stream.write("     values, and `unplaced_id` key with string value.\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this usage message\n")
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  1. The in.tsv file is an Orthogroups.tsv file with header defined.\n")
    stream.write("\n")
    stream.write("  2. ClusterVenn3 requires `.txt` file suffix\n")
    stream.write("\n")
    stream.write("  3. Flag order determines execution order for the following optional flags:\n")
    stream.write("     `--output-reference-names`, `--input-reference-names`,\n")
    stream.write("     `--prepend-sample-prefix`, `--remove-sample-prefix`,\n")
    stream.write("     `--samples-file`, `--samples-list`, and `--sort-by`\n")
    stream.write("\n\n%s" % message)
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    sys.exit(exitcode)
    

    
def main(argv):
    short_flags = 'D:d:hIiNO:o:PrS:s:uvy:'
    long_flags = (
        'replace-delim',
        'prefix-delim=',
        'ignore-unplaced-strictly',
        'ignore-unplaced-leniently',
        'output-reference-names',
        'output-file=', 'outfile=',
        'output-type=', 'orthovenn',
        'samples-list=','sample-order=','species-order=',
        'samples-file=','sample-order-file=','species-order-file=',
        'prepend-sample-prefix','prefix-sample-names','prefix-species-names',
        'remove-sample-prefix','remove-sample-names','remove-species-names',
        'ignore-unlocalized',
        'sort-by=',
        'yaml=',
        'help'
    )
    try:
        options, arguments = getopt(argv, short_flags, long_flags)
    except GetoptError as message:
        usage(message)
                 
    output_type = ENUM_OUTPUT.ORTHOFINDER
    config_file = None
    output_file = STDIO
    prefix_delim = '|'
    replace_delim = prefix_delim
    ignore_unplaced = 0
    ignore_unlocalized = _placed
    input_reference_names = False
    output_reference_names = False
    call_order = []
    ogid_as_sample = False
    for flag, value in options:
        if  flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-D','--replace-delim'):
            replace_delim = value        
        elif flag in ('-d','--prefix-delim'):
            prefix_delim = value
        elif flag in ('-I','ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-i','ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-N','--output-reference-names'):
            call_order.append((ENUM_CALLKEY.OUTPUT_REFS, None))
            call_order.append((ENUM_CALLKEY.INPUT_REFS, None))
            output_reference_names = True
        elif flag in ('-O','--output-type'):
            output_type = value
            if output_type not in ENUM_OUTPUT:
                usage("Invalid enumerative `%s` value: `%s`" % (flag, value))
        elif flag in ('-o','--outfile','--output-file'):
            output_file = value
        elif flag in ('-P','--remove-sample-prefix',
                      '--remove-sample-names','--remove-species-names'):
            call_order.append((ENUM_CALLKEY.REMOVE_SAMPLES, None))
        elif flag in ('-p','--prepend-sample-prefix',
                      '--prefix-sample-names','--prefix-species-names'):
            call_order.append((ENUM_CALLKEY.PREFIX_SAMPLES, None))
        elif flag in ('-r','--input-reference-names'):
            call_order.append((ENUM_CALLKEY.INPUT_REFS, None))
            input_reference_names = True
        elif flag in ('-S','--samples-file','--sample-order-file','--species-order-file'):
            sample_names = read_listfiles(value)
            has_ogid_as_sample = has_ogid(sample_names)
            if has_ogid_as_sample:
                ogid_as_sample = has_ogid_as_sample
            call_order.append((ENUM_CALLKEY.SUBSET_SAMPLES, sample_names))
        elif flag in ('-s','--samples-list','--sample-order','--species-order'):
            sample_names = value.strip(_COMMA).split(_COMMA)
            has_ogid_as_sample = has_ogid(sample_names)
            if has_ogid_as_sample:
                ogid_as_sample = has_ogid_as_sample
            call_order.append((ENUM_CALLKEY.SUBSET_SAMPLES, sample_names))
        elif flag in ('--sort-by',):
            sortby = int(value)
            if not (abs(sortby) & (ENUM_SORTBY.OGID|ENUM_SORTBY.MEMBERS|ENUM_SORTBY.SAMPLES)):
                 usage("Invalid enumerative `%s` value: `%s`" % (flag, value))
            call_order.append((ENUM_CALLKEY.SORTBY, sortby))
        elif flag in ('-u','--ignore-unlocalized'):
            ignore_unlocalized = _localized
        elif flag in ('-v','--orthovenn'):
            output_type = ENUM_OUTPUT.ORTHOVENN
        elif flag in ('-y','--yaml'):
            config_file = value

    if input_reference_names or output_reference_names:
        if config_file is None:
            usage('`--yaml` config file required')
            
    if num(arguments) != 1:
        usage("Unexpected number of arguments")

    view(
        arguments[0],
        config_file=config_file,
        call_order=call_order,
        ogid_as_sample=ogid_as_sample,
        prefix_delim=prefix_delim,
        replace_delim=replace_delim,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        output_reference_names=output_reference_names,
        output_type=output_type,
        output_file=output_file
    )
                 
    return 0



def view(
        ortho_file,
        config_file=None,
        call_order=(),
        ogid_as_sample=False,
        prefix_delim='|',
        replace_delim='|',
        ignore_unplaced=None,
        ignore_unlocalized=False,
        input_reference_names=False,
        output_reference_names=False,
        output_type=ENUM_OUTPUT.ORTHOFINDER,
        output_file=None
):
    ortho = OrthogroupsFile(ortho_file)
    config = None
    if config_file:
        config = SampleConfigFile(
            config_file,
            load_files=True,
            map_assigned_molecule=True
        )
                 
    ortho = api.view(
        ortho, config=config,
        call_order=call_order,
        prefix_delim=prefix_delim,
        replace_delim=replace_delim,
        ignore_unplaced=ignore_unplaced,
        ignore_unlocalized=ignore_unlocalized,
        input_reference_names=input_reference_names,
        output_reference_names=output_reference_names
    )
    if output_type == ENUM_OUTPUT.ORTHOVENN:
        with open(output_file, 'w') as output:
            for record in OrthoVennFormatter(ortho):
                output.write(record + _EOL)
            
    elif output_type == ENUM_OUTPUT.ANCHORS:
        if ogid_as_sample:
            if num(ortho.samples) != 1:
                raise Exception(_TWO_SAMPLES_REQUIRED)
            sample1 = ogid_as_sample
            sample2 = ortho.samples[0]
        elif num(ortho.samples) == 2:
            sample1 = ortho.samples[0]
            sample2 = ortho.samples[1]
        else:
            raise Exception(_TWO_SAMPLES_REQUIRED)

        with open(output_file, 'w') as output:
            for record in MCScanAnchorsFormatter(ortho, sample1, sample2):
                output.write(record + _EOL)
            
    elif output_type == ENUM_OUTPUT.ORTHOFINDER:
        ortho.to_file(output_file)
    
    else:
        raise AssertionError("Invalid `--output-type`")
