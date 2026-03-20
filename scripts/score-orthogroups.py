#!/usr/bin/env python


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
    for group_index in range(num(ortho.groups)):
        for species_index in range(num(ortho.species)):
            if ortho.groups[group_index][species_index] is None:
                continue
            for member in ortho.groups[group_index][species_index]:
                if member not in joint_freq:
                    joint_freq[member] = [0] * num_clusters
                    
                joint_freq[member][cluster_index[ortho.clusters[group_index]]] += ortho.counts[group_index]
                total += ortho.counts[group_index]

    total = float(total)
    for member in joint_freq:
        for clust_index in range(num(joint_freq[member])):
            joint_freq[member][clust_index] = float(joint_freq[member][clust_index]) / total

    return joint_freq



def _calc_marginal_prob(ortho):
    total = 0
    mrgnl_freq = {}
    for group_index in range(num(ortho.groups)):
        for species_index in range(num(ortho.species)):
            if ortho.groups[group_index][species_index] is None:
                continue
            for member in ortho.groups[group_index][species_index]:
                if member in mrgnl_freq:
                    mrgnl_freq[member] += ortho.counts[group_index]
                else:
                    mrgnl_freq[member] = ortho.counts[group_index]
                total += ortho.counts[group_index]

    total = float(total)
    for member in mrgnl_freq:
        mrgnl_freq[member] = float(mrgnl_freq[member]) / total
    
    return mrgnl_freq



def calc_cluster_joint_freqs(
        ortho, namemap, ignore_groups=set(),
        ignore_unplaced=0, is_placed=_placed, init=1e-6
):
    cluster_freqs = dict()
    cluster_names = sorted(set(ortho.clusters))
    cluster_inits = {species_name: {None: 1} for species_name in ortho.species}
    for species_name in ortho.species:
        if species_name not in namemap.species:
            raise KeyError("Species name not found in config: %s" % species_name)
        for member in namemap.species[species_name].references:
            if member == namemap.species[species_name].unplaced_id or \
               not is_placed(namemap.species[species_name].references[member]):
                if ignore_unplaced == _STRICT:
                    continue
            cluster_inits[species_name][member] = 1

    for cluster_name in cluster_names:
        if cluster_name in ignore_groups:
            continue
        cluster_freqs[cluster_name] = {
            species_name: cluster.copy() for species_name, cluster in cluster_inits.items()
        }
        
    for group_index in range(num(ortho.groups)):
        if ortho.clusters[group_index] in ignore_groups:
            continue
        for species_index in range(num(ortho.species)):
            n = 0
            for member in ortho.groups[group_index][species_index]:
                if member == namemap.species[ortho.species[species_index]].unplaced_id or \
                   not is_placed(namemap.species[ortho.species[species_index]].references[member]):
                    if ignore_unplaced == _STRICT:
                        continue
                cluster_freqs[ortho.clusters[group_index]][ortho.species[species_index]][member] += ortho.counts[group_index]
                n += 1
            if n == 0:
                cluster_freqs[ortho.clusters[group_index]][ortho.species[species_index]][None] += 1

    
    for cluster_name in cluster_names:
        if cluster_name in ignore_groups:
            continue        
        for species_name in cluster_freqs[cluster_name]:
            total = float(sum(cluster_freqs[cluster_name][species_name].values()))

            cluster_freqs[cluster_name][species_name] = {
                key: val / total for key, val in cluster_freqs[cluster_name][species_name].items()
            }
            cluster_freqs[cluster_name][species_name][None] = init
            
    return cluster_freqs



def calc_cluster_marginal_freqs(ortho, ignore_groups=set()):
    total = 0
    cluster_freqs = dict()
    cluster_names = sorted(set(ortho.clusters))
    for group_index in range(num(ortho.groups)):
        if ortho.clusters[group_index] in ignore_groups:
            continue
        try:
            cluster_freqs[ortho.clusters[group_index]] += ortho.counts[group_index]
        except KeyError:
            cluster_freqs[ortho.clusters[group_index]]  = ortho.counts[group_index]
        total += ortho.counts[group_index]
        
    cluster_freqs = {key: val / total for key, val in cluster_freqs.items()}
        
    return cluster_freqs



def calc_conditional_prob(
        orthoA, orthoB, namemap, ignore_groups=set(),
        ignore_unplaced=0, is_placed=_placed, init=1e-6
):
    """
    Calculate the joint probability of cluster (A) and chr (B), divided by the
    marginal probability of a chr.

    By Bayes' Theorum:

      P(A | B) =  (P(B | A) * P(A)) / P(B)
               = ((P(B n A) / P(A)) * P(A)) / P(B)
               =   P(B n A) / P(B)
    """
    joint_prob = \
        calc_cluster_joint_freqs(
            orthoB,
            namemap,
            ignore_groups=ignore_groups,
            ignore_unplaced=ignore_unplaced,
            is_placed=_placed
        )  #_calc_joint_prob(orthoB)
    
    marginal_prob = \
        calc_cluster_marginal_freqs(
            orthoB,
            ignore_groups=ignore_groups,
        )  # _calc_marginal_prob(orthoB)

    probs = ProbabilitiesTable()
    probs.colnames = sorted(set(orthoB.clusters))
    probs.colindex = index_list(probs.colnames)
    prior_species = set(orthoB.species)
    num_clusters = num(probs.colnames)
    
    for group_index in range(num(orthoA.groups)):
        cluster_prob = [-1] * num_clusters

        for cluster_name in joint_prob:
            if cluster_name in ignore_groups:
                continue
            
            for species_index in range(num(orthoA.species)):
                if orthoA.species[species_index] not in joint_prob[cluster_name]:
                    continue
                n = 0
                prob = 0.0
                for member in orthoA.groups[group_index][species_index]:
                    if member == namemap.species[orthoA.species[species_index]].unplaced_id or \
                       not is_placed(namemap.species[orthoA.species[species_index]].references[member]):
                        if ignore_unplaced == _STRICT:
                            continue
                    prob += joint_prob[cluster_name][orthoA.species[species_index]][member]
                    n += 1
                if n == 0:
                    prob = joint_prob[cluster_name][orthoA.species[species_index]][None]
                    
                if cluster_prob[probs.colindex[cluster_name]] < 0:
                    cluster_prob[probs.colindex[cluster_name]]  = prob
                else:
                    cluster_prob[probs.colindex[cluster_name]] *= prob
                    
            cluster_prob[probs.colindex[cluster_name]] *= marginal_prob[cluster_name]
            
        cluster_prob = list(map(_min0, cluster_prob))
        cluster_norm = sum(cluster_prob) or 1.0
        cluster_prob = list(map(lambda prob: prob / cluster_norm, cluster_prob))
        
        probs.rownames.append(orthoA.ids[group_index])
        probs.matrix.append(cluster_prob)

    probs.rowindex = index_list(probs.rownames)
        
    return probs



def assign_clusters(ortho, probs):
    cluster_indices = list(range(num(probs.colnames)))
    for group_index in range(num(ortho.groups)):
        group_id = ortho.ids[group_index]
        if group_id in probs.rowindex:
            cluster_id = max(
                cluster_indices,
                key=probs.matrix[probs.rowindex[group_id]].__getitem__
            )
            ortho.clusters[group_index] = probs.colnames[cluster_id]
            ortho.probabilities[group_index] = probs.matrix[probs.rowindex[group_id]][cluster_id]

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
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80    
    stream.write("  -E,--check-errors\n")
    stream.write("     Check for classification errors in classified.tsv\n")
    stream.write("\n")
    stream.write("  -G,--ignore-groups <name1>[,<name2>[,...]]\n")
    stream.write("     Exclude groups, designated by their comma-separated list of group\n")
    stream.write("     names, from forming Bayesian classification groups. Orthogroups in\n")
    stream.write("     the designated synteny groups will be clustered to other groups.\n")
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
    short_options = 'hEG:iIuP:'
    long_options = (
        'help',
        'check-errors',
        'ignore-groups=',
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
    ignore_groups = set()
    ignore_unplaced = 0  # _STRICT
    output_posterior_file = None
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-E','--check-errors'):
            check_errors = True
        elif flag in ('-u','--ignore-unlocalized'):
            is_placed = _localized
        elif flag in ('-i','--ignore-unplaced-leniently'):
            ignore_unplaced = _LENIENT
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-P','--output-posterior-prob-file'):
            output_posterior_file = value
        elif flag in ('-G','--ignore-groups'):
            ignore_groups = set(value.strip(_COMMA).split(_COMMA))

            
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
    
    pprobs = calc_conditional_prob(orthoU, orthoM, config, ignore_groups, ignore_unplaced, is_placed)

    orthoU = assign_clusters(orthoU, pprobs)

    if check_errors:
        col_index = pprobs.colindex
        row_index = pprobs.rowindex
        errors = ClusterErrorOrthogroups()
        errors.species = orthoU.species
        for group_index in range(num(orthoU.groups)):
            if orthoU.clusters[group_index] != orthoM.clusters[group_index]:
                errors.clusters.append((
                    orthoM.clusters[group_index],
                    orthoU.clusters[group_index]
                ))
                errors.probabilities.append((
                    pprobs.matrix[row_index[orthoM.ids[group_index]]][col_index[orthoM.clusters[group_index]]],
                    pprobs.matrix[row_index[orthoU.ids[group_index]]][col_index[orthoU.clusters[group_index]]]
                ))
                errors.counts.append(orthoU.counts[group_index])
                errors.groups.append(orthoU.groups[group_index])
                errors.ids.append(orthoU.ids[group_index])

        orthoU = errors

    if output_posterior_file is not None:
        pprobs.to_table(output_posterior_file)
    orthoU.to_table()



if __name__ == '__main__':
    main(sys.argv[1:])
