#!/usr/bin/env python

# TODO: modify to allow reading different file input Orthofinder-like file formats

import os

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Use manual clustering to classify new groups'


import io
import sys
import getopt

from og.core.members import _LENIENT, _STRICT
from og.core.compression import is_stream, open
from og.core.parsers.config import SpeciesConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.orthogroups import CountedClusteredOrthogroups
from og.core.parsers.orthogroups import ClusterErrorOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
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


class ProbabilitiesTable(object):
    def __init__(self):
        self.colnames = []
        self.rownames = []
        self.colindex = {}
        self.rowindex = {}
        self.matrix = []

    def format_header(self):
        return "Cluster\t%s" % _TAB.join(self.colnames)

    def format_record(self, index=None):
        if index is None:
            raise ValueError("no index given")
        return "%s\t%s" % (
            self.rownames[index],
            _TAB.join(map("{:.3e}".format, self.matrix[index]))
        )
        
    def to_table(self, infile=sys.stdout, **kwargs):
        if is_stream(infile):
            stream = infile
            close = False
        else:
            if 'mode' in kwargs:
                if 'r' in kwargs['mode']:
                    raise ValueError("%s.to_table() is write-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'w'
                stream = open(infile, **kwargs)
                close = True

        stream.write(self.format_header() + _EOL)
        for i in range(num(self.rownames)):
            stream.write(self.format_record(index=i) + _EOL)
        if close:
            stream.close()



def _min0(x):
    return 0.0 if x < 0.0 else x


def index_list(lst):
    return { item: i for i, item in enumerate(lst) }

    
def _calc_joint_prob(ortho):
    total = 0
    joint_freq = {}
    num_clusters = num(set(ortho.clusters))
    cluster_index = index_list(sorted(set(ortho.clusters)))
    for g in range(num(ortho.groups)):
        for s in range(num(ortho.species)):
            if ortho.groups[g][s] is None:
                continue
            for member in ortho.groups[g][s]:
                if member not in joint_freq:
                    joint_freq[member] = [0] * num_clusters
                    
                joint_freq[member][cluster_index[ortho.clusters[g]]] += ortho.counts[g]
                total += ortho.counts[g]

    total = float(total)
    for member in joint_freq:
        for c in range(num(joint_freq[member])):
            joint_freq[member][c] = float(joint_freq[member][c]) / total

    return joint_freq


def _calc_marginal_prob(ortho):
    total = 0
    mrgnl_freq = {}
    for g in range(num(ortho.groups)):
        for s in range(num(ortho.species)):
            if ortho.groups[g][s] is None:
                continue
            for member in ortho.groups[g][s]:
                if member in mrgnl_freq:
                    mrgnl_freq[member] += ortho.counts[g]
                else:
                    mrgnl_freq[member] = ortho.counts[g]
                total += ortho.counts[g]

    total = float(total)
    for member in mrgnl_freq:
        mrgnl_freq[member] = float(mrgnl_freq[member]) / total
    
    return mrgnl_freq

    
def calc_conditional_prob(orthoA, orthoB, namemap, is_placed=_placed, ignore_unplaced=0, init=1e-6):
    """
    Calculate the joint probability of chr and cluster, divided by the marginal
    probability of a chr.

    By Bayes' Theorum:

      P(A | B) =  (P(B | A) * P(A)) / P(B)
               = ((P(B n A) / P(A)) * P(A)) / P(B)
               =   P(B n A) / P(B)
    """
    joint_prob = _calc_joint_prob(orthoB)
    marginal_prob = _calc_marginal_prob(orthoB)
    
    probs = ProbabilitiesTable()
    probs.colnames = sorted(set(orthoB.clusters))
    prior_species = set(orthoB.species)    
    num_clusters = num(probs.colnames)
    for g in range(num(orthoA.groups)):
        cluster_prob = [-1] * num_clusters
        
        for s in range(num(orthoA.species)):
            if orthoA.species[s] not in prior_species:
                continue
            for c in range(num_clusters):
                prob = init
                count = 0
                for member in orthoA.groups[g][s]:
                    if member == namemap.species[orthoA.species[s]].unplaced_id or \
                       not is_placed(namemap.species[orthoA.species[s]].references[member]):
                        if ignore_unplaced == _STRICT:
                            continue
                    try:
                        prob += joint_prob[member][c] / marginal_prob[member]
                        count += 1
                    except KeyError:
                        pass
                    
                if cluster_prob[c] < 0:
                    cluster_prob[c]  = (prob / float(max(1, count)))
                else:
                    cluster_prob[c] *= (prob / float(max(1, count)))

        cluster_prob = list(map(_min0, cluster_prob))
        cluster_norm = sum(cluster_prob) or 1.0
        cluster_prob = list(map(lambda p: p/cluster_norm, cluster_prob))
        
        probs.rownames.append(orthoA.ids[g])
        probs.matrix.append(cluster_prob)

    probs.colindex = index_list(probs.colnames)
    probs.rowindex = index_list(probs.rownames)
        
    return probs


def assign_clusters(ortho, probs):
    cindices = list(range(num(probs.colnames)))
    for g in range(num(ortho.groups)):
        gid = ortho.ids[g]
        if gid in probs.rowindex:
            cid = max(
                cindices,
                key=probs.matrix[probs.rowindex[gid]].__getitem__
            )
            ortho.clusters[g] = probs.colnames[cid]
            ortho.probabilities[g] = probs.matrix[probs.rowindex[gid]][cid]

    return ortho


def _copy_to_CountedClusteredOrthogroups(ortho):
    if isinstance(ortho, CountedClusteredOrthogroups):
        return ortho
    clust = CountedClusteredOrthogroups()
    clust.species = ortho.species
    clust.ids = ortho.ids
    clust.clusters = ortho.ids.copy()
    clust.groups = ortho.groups
    clust.counts = [1] * num(ortho.groups)
    clust.probabilities = [-1] * num(ortho.groups)
    return clust


def open_inferred_format(filename):
    filehandle = open(filename, 'rt')

    seekable = False

    firstline = next(filehandle)
    if firstline.startswith('Count\tCluster'):
        constructor = CountedClusteredOrthogroups
    else:
        constructor = OrthoFinderOrthogroups

    if filehandle.seekable():
        filehandle.seek(0)
        return _copy_to_CountedClusteredOrthogroups(constructor(filehandle))
    else:
        ortho = constructor()
        ortho.from_string(firstline + ''.join(filehandle))
        return _copy_to_CountedClusteredOrthogroups(ortho)

    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <classified.tsv> [unclassified.tsv] <in.yaml>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -E,--check-errors\n")
    stream.write("     Check for classification errors in classified.tsv\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     Strictly ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains no chromosomal sequences, that cell\n")
    stream.write("     then contains no members.\n")
    stream.write("\n")
    stream.write("  -i,--ignore-unplaced-leniently\n")
    stream.write("     Leniently ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains only unplaced (ie, non-chomosomal)\n")
    stream.write("     sequences, that cell contains members.\n")
    stream.write("\n")
    stream.write("  -P,--output-posterior-prob-file <file>\n")
    stream.write("     Write the Bayesian posterior probabilities table to file.\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized sequences to their\n")
    stream.write("     designated sequence names, not to their placed chromosome names.\n")
    stream.write("\n")    
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  classified.tsv is a mrg.counts.tsv file output by cluster-orthogroups.py\n")
    stream.write("  or an OrthoFinder orthogroups.tsv-formatted files. Lines belonging to the\n")
    stream.write("  same syntenic groups are expected to be grouped together and such groups\n")
    stream.write("  preceded by `##group=N` meta lines to demarcate the group and define their\n")
    stream.write("  cluster ID (N).\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    
def main(argv):
    short_options = 'hEiIuP:'
    long_options = (
        'help',
        'check-errors',
        'ignore-unlocalized',
        'ignore-unplaced-strictly',
        'ignore-unplaced-leniently',
        'output-posterior-prob-file='
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    check_errors = False
    is_placed = _placed
    ignore_unplaced = 0  # _STRICT
    output_posterior_file = None
    for flag, value in options:
        if   flag in ('-h','--help'): usage(exitcode=0)
        elif flag in ('-E','--check-errors'): check_errors = True
        elif flag in ('-u','--ignore-unlocalized'): is_placed = _localized
        elif flag in ('-i','--ignore-unplaced-leniently'): ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'): ignore_unplaced = _STRICT
        elif flag in ('-P','--output-posterior-prob-file'): output_posterior_file = value

    if num(arguments) != 2 and \
       num(arguments) != 3:
        usage('Unexpected number of arguments')

    if num(arguments) == 2:
        arguments.append(arguments[1])
        arguments[1] = arguments[0]
    if check_errors:
        arguments[1] = arguments[0]

        
    orthoM = CountedClusteredOrthogroups(arguments[0])
    orthoU = open_inferred_format(arguments[1])
    config = SpeciesConfigFile(arguments[2], load_files=True, map_assigned_molecule=True)
    
    pprobs = calc_conditional_prob(orthoU, orthoM, config, is_placed, ignore_unplaced)

    orthoU = assign_clusters(orthoU, pprobs)

    if check_errors:
        cindex = pprobs.colindex
        rindex = pprobs.rowindex
        errors = ClusterErrorOrthogroups()
        errors.species = orthoU.species
        for g in range(num(orthoU.groups)):
            if orthoU.clusters[g] != orthoM.clusters[g]:
                errors.clusters.append((
                    orthoM.clusters[g],
                    orthoU.clusters[g]
                ))
                errors.probabilities.append((
                    pprobs.matrix[rindex[orthoM.ids[g]]][cindex[orthoM.clusters[g]]],
                    pprobs.matrix[rindex[orthoU.ids[g]]][cindex[orthoU.clusters[g]]]
                ))
                errors.counts.append(orthoU.counts[g])
                errors.groups.append(orthoU.groups[g])
                errors.ids.append(orthoU.ids[g])

        orthoU = errors

    if output_posterior_file is not None:
        pprobs.to_table(output_posterior_file)
    orthoU.to_table()

    
if __name__ == '__main__':
    main(sys.argv[1:])
