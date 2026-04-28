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

from og.core.utils import index_list
from og.core.members import _LENIENT, _STRICT
from og.core.compression import is_stream, open
from og.core.parsers.config import SampleConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.orthogroups import CountedClusteredOrthogroups
from og.core.parsers.orthogroups import ClusterErrorOrthogroups
from og.core.parsers.orthogroups import _COMPRESSION_FLAGS
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
        
    def to_file(self, file=sys.stdout, **kwargs):
        if is_stream(file):
            stream = file
            close = False
        else:
            if 'mode' in kwargs:
                if 'r' in kwargs['mode']:
                    raise ValueError("%s.to_file() is write-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'w'
            kwargs = {
                k:v for k,v in kwargs.items() if k in _COMPRESSION_FLAGS
            }
            stream = open(file, **kwargs)
            close = True

        stream.write(self.format_header() + _EOL)
        for i in range(num(self.rownames)):
            stream.write(self.format_record(index=i) + _EOL)
        if close:
            stream.close()

    to_table = to_file
    

    
def _min0(x):
    return 0.0 if x < 0.0 else x


def _calc_joint_prob(ortho):
    total = 0
    joint_freq = {}
    cluster_index = index_list(sorted(set(g.cluster for g in ortho.groups)))
    num_clusters = num(cluster_index)
    for group in ortho.groups:
        for sample in ortho.samples:
            if not group[sample.index]:
                continue
            for member in group[sample.index]:
                if member not in joint_freq:
                    joint_freq[member] = [0] * num_clusters
                    
                joint_freq[member][cluster_index[group.cluster]] += group.count
                total += group.count

    total = float(total)
    for member in joint_freq:
        for clust_index in range(num(joint_freq[member])):
            joint_freq[member][clust_index] = float(joint_freq[member][clust_index]) / total

    return joint_freq



def _calc_marginal_prob(ortho):
    total = 0
    mrgnl_freq = {}
    for group in ortho.groups:
        for sample in ortho.samples:
            if not group[sample.index]:
                continue
            for member in group[sample.index]:
                if member in mrgnl_freq:
                    mrgnl_freq[member] += group.count
                else:
                    mrgnl_freq[member] = group.count
                total += group.count

    total = float(total)
    for member in mrgnl_freq:
        mrgnl_freq[member] = float(mrgnl_freq[member]) / total
    
    return mrgnl_freq



def calc_cluster_joint_freqs(
        ortho, namemap, ignore_groups=set(),
        ignore_unplaced=0, is_placed=_placed, init=1e-6
):
    cluster_freqs = dict()
    cluster_names = sorted(set(g.cluster for g in ortho.groups))
    cluster_inits = {s.id: {None: 1} for s in ortho.samples}
    for sample in ortho.samples:
        if sample.id not in namemap.samples:
            raise KeyError(
                "Sample not found in YAML file: '%s'" % str(sample.id)
            )
        for member in namemap.samples[sample.id].references:
            if member == namemap.samples[sample.id].unplaced_id or \
               not is_placed(namemap.samples[sample.id].references[member]):
                if ignore_unplaced == _STRICT:
                    continue
            cluster_inits[sample.id][member] = 1

    for cluster_name in cluster_names:
        if cluster_name in ignore_groups:
            continue
        cluster_freqs[cluster_name] = {
            sample_id: cluster.copy() for sample_id, cluster in cluster_inits.items()
        }
        
    for group in ortho.groups:
        if group.cluster in ignore_groups:
            continue
        for sample in ortho.samples:
            n = 0
            for member in group[sample.index]:
                if member == namemap.samples[sample.id].unplaced_id or \
                   not is_placed(namemap.samples[sample.id].references[member]):
                    if ignore_unplaced == _STRICT:
                        continue
                if member not in cluster_freqs[group.cluster][sample.id]:
                    continue
                
                cluster_freqs[group.cluster][sample.id][member] += group.count
                n += 1
                
            if n == 0:
                cluster_freqs[group.cluster][sample.id][None] += 1

    
    for cluster_name in cluster_names:
        if cluster_name in ignore_groups:
            continue        
        for sample_name in cluster_freqs[cluster_name]:
            total = float(sum(cluster_freqs[cluster_name][sample_name].values()))

            cluster_freqs[cluster_name][sample_name] = {
                key: val / total for key, val in cluster_freqs[cluster_name][sample_name].items()
            }
            cluster_freqs[cluster_name][sample_name][None] = init
            
    return cluster_freqs



def calc_cluster_marginal_freqs(ortho, ignore_groups=set()):
    total = 0
    cluster_freqs = dict()
    cluster_names = sorted(set(g.cluster for g in ortho.groups))
    for group in ortho.groups:
        if group.cluster in ignore_groups:
            continue
        try:
            cluster_freqs[group.cluster] += group.count
        except KeyError:
            cluster_freqs[group.cluster] = group.count
        total += group.count
        
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
    probs.colnames = sorted(set(g.cluster for g in orthoB.groups)) + [None]
    probs.colindex = index_list(probs.colnames)
    prior_samples = set(orthoB.samples)
    num_clusters = num(probs.colnames)
    
    for group in orthoA.groups:
        cluster_prob = [-1] * num_clusters

        for cluster_name in joint_prob:
            if cluster_name in ignore_groups:
                continue
            
            for sample in orthoA.samples:
                if sample.id not in joint_prob[cluster_name]:
                    continue
                if not group[sample.index]:
                    continue
                n = 0
                prob = 0.0
                for member in group[sample.index]:
                    if member not in joint_prob[cluster_name][sample.id]:
                        continue                    
                    if member == namemap.samples[sample.id].unplaced_id or \
                       not is_placed(namemap.samples[sample.id].references[member]):
                        if ignore_unplaced == _STRICT:
                            continue
                    prob += joint_prob[cluster_name][sample.id][member]
                    n += 1
                if n == 0:
                    prob = joint_prob[cluster_name][sample.id][None]
                    
                if cluster_prob[probs.colindex[cluster_name]] < 0:
                    cluster_prob[probs.colindex[cluster_name]]  = prob
                else:
                    cluster_prob[probs.colindex[cluster_name]] *= prob
                    
            cluster_prob[probs.colindex[cluster_name]] *= marginal_prob[cluster_name]
            
        cluster_prob = list(map(_min0, cluster_prob))
        cluster_norm = sum(cluster_prob) or 1.0
        cluster_prob = list(map(lambda prob: prob / cluster_norm, cluster_prob))
        
        probs.rownames.append(group.id)
        probs.matrix.append(cluster_prob)

    probs.rowindex = index_list(probs.rownames)
        
    return probs



def assign_clusters(ortho, probs):
    cluster_indices = list(range(num(probs.colnames)))
    for group in ortho.groups:
        if group.id in probs.rowindex:
            cluster_id = max(
                cluster_indices,
                key=probs.matrix[probs.rowindex[group.id]].__getitem__,
                default=None
            )
            group.cluster = probs.colnames[cluster_id]
            group.probability = probs.matrix[probs.rowindex[group.id]][cluster_id]

    return ortho



def _copy_to_CountedClusteredOrthogroups(old_ortho):
    if isinstance(old_ortho, CountedClusteredOrthogroups):
        return old_ortho
    
    new_ortho = CountedClusteredOrthogroups(samples=old_ortho.samples)
    for old_group in old_ortho.groups:
        new_group = new_ortho.new_group(append=True)
        for i in range(num(old_ortho.samples)):
            new_group[i] = old_group[i]
        new_group.id = old_group.id
        new_group.cluster = old_group.id
        new_group.count = 1
        new_group.probability = -1
    return new_ortho



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
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
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
    short_options = 'hEG:iIuP:o:'
    long_options = (
        'help',
        'check-errors',
        'ignore-groups=',
        'ignore-unlocalized',
        'ignore-unplaced-strictly',
        'ignore-unplaced-leniently',
        'output-posterior-prob-file=',
        'output-file=',
    )
    try:
        options, arguments = getopt.getopt(argv, short_options, long_options)
    except getopt.GetoptError as message:
        usage(message)

    check_errors = False
    is_placed = _placed
    output_file = sys.stdout
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
        elif flag in ('-o','--output-file'):
            output_file = open(value, 'w')
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
    config = SampleConfigFile(arguments[2], load_files=True, map_assigned_molecule=True)
    
    pprobs = calc_conditional_prob(orthoU, orthoM, config, ignore_groups, ignore_unplaced, is_placed)

    orthoU = assign_clusters(orthoU, pprobs)

    if check_errors:
        col_index = pprobs.colindex
        row_index = pprobs.rowindex
        errors = ClusterErrorOrthogroups(samples=orthoU.samples)
        for group_index in range(num(orthoU.groups)):
            if orthoU.groups[group_index].cluster != orthoM.groups[group_index].cluster:
                error = errors.new_group(append=True)
                for i in range(num(orthoU.samples)):
                    error[i] = orthoU.groups[group_index][i]
                error.id = orthoU.groups[group_index].id
                error.count = orthoU.groups[group_index].count
                error.cluster = (
                    orthoM.groups[group_index].cluster,
                    orthoU.groups[group_index].cluster
                )
                error.probability = (
                    pprobs.matrix[row_index[orthoM.groups[group_index].id]][col_index[orthoM.groups[group_index].cluster]],
                    pprobs.matrix[row_index[orthoU.groups[group_index].id]][col_index[orthoU.groups[group_index].cluster]]
                )

        orthoU = errors

    if output_posterior_file is not None:
        pprobs.to_file(output_posterior_file)
    orthoU.to_file(output_file)



if __name__ == '__main__':
    main(sys.argv[1:])
