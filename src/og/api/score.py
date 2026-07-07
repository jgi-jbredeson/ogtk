
import sys

from og.core.utils import index_list
from og.core.members import _STRICT
from og.core.compressio import is_stream, open
from og.core.assembly_report import is_placed as _placed
from og.core.orthogroups.parsers import ClusterErrorOrthogroups
from og.core.orthogroups.parsers import _COMPRESSION_FLAGS
from og.constants import _TAB, range


num = len


class ProbabilitiesTable(object):
    def __init__(self):
        self.colnames = []
        self.rownames = []
        self.colindex = {}
        self.rowindex = {}
        self.matrix = []

    def format_header(self):
        return "Cluster\t%s" % _TAB.join(map(str, self.colnames))

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

        print(self.format_header(), file=stream)
        for i in range(num(self.rownames)):
            print(self.format_record(index=i), file=stream)
        if close:
            stream.close()

    to_table = to_file
    

    
def _min0(x):
    return 0.0 if x < 0.0 else x



def _calc_joint_probs(ortho):
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



def _calc_marginal_probs(ortho):
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
        ortho, config, ignore_groups=set(),
        ignore_unplaced=0, ignore_unlocalized=_placed, init=1e-6
):
    config.check_samples([s.id for s in ortho.samples])
    
    cluster_freqs = dict()
    cluster_names = sorted(set(g.cluster for g in ortho.groups))
    cluster_inits = {s.id: {None: 1} for s in ortho.samples}
    for sample in ortho.samples:
        for member in config.samples[sample.id].references:
            if member == config.samples[sample.id].unplaced_id or \
               not ignore_unlocalized(config.samples[sample.id].references[member]):
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
                if member == config.samples[sample.id].unplaced_id or \
                   not ignore_unlocalized(config.samples[sample.id].references[member]):
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



def calc_conditional_probs(
        classify_ortho, training_ortho, config, ignore_groups=set(),
        ignore_unplaced=0, ignore_unlocalized=_placed, init=1e-6
):
    """
    Calculate the joint probability of cluster (A) and chr (B), divided by the
    marginal probability of a chr.

    By Bayes' Theorum:

      P(A | B) =  (P(B | A) * P(A)) / P(B)
               = ((P(B n A) / P(A)) * P(A)) / P(B)
               =   P(B n A) / P(B)
    """
    config.check_samples([s.id for s in classify_ortho.samples])
    config.check_samples([s.id for s in training_ortho.samples])
    
    joint_prob = \
        calc_cluster_joint_freqs(
            training_ortho,
            config,
            ignore_groups=ignore_groups,
            ignore_unplaced=ignore_unplaced,
            ignore_unlocalized=ignore_unlocalized
        )  #_calc_joint_prob(training_ortho)
    
    marginal_prob = \
        calc_cluster_marginal_freqs(
            training_ortho,
            ignore_groups=ignore_groups,
        )  # _calc_marginal_prob(training_ortho)

    probs = ProbabilitiesTable()
    probs.colnames = sorted(set(g.cluster for g in training_ortho.groups)) + [None]
    probs.colindex = index_list(probs.colnames)
    num_clusters = num(probs.colnames)
    
    for group in classify_ortho.groups:
        cluster_prob = [-1] * num_clusters

        for cluster_name in joint_prob:
            if cluster_name in ignore_groups:
                continue
            
            for sample in classify_ortho.samples:
                if sample.id not in joint_prob[cluster_name]:
                    continue
                if not group[sample.index]:
                    continue
                n = 0
                prob = 0.0
                for member in group[sample.index]:
                    if member not in joint_prob[cluster_name][sample.id]:
                        continue                    
                    if member == config.samples[sample.id].unplaced_id or \
                       not ignore_unlocalized(config.samples[sample.id].references[member]):
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
            cluster_index = max(
                cluster_indices,
                key=probs.matrix[probs.rowindex[group.id]].__getitem__,
                default=None
            )
            group.cluster = probs.colnames[cluster_index]
            group.probability = probs.matrix[probs.rowindex[group.id]][cluster_index]

    return ortho


    
def score(    # orthoU  orthoM
        classify_ortho, training_ortho, config,
        check_errors=False, ignore_groups=(),
        ignore_unplaced=False, ignore_unlocalized=_placed
):
    pprobs = \
        calc_conditional_probs(
            classify_ortho,  # orthoU
            training_ortho,  # orthoM
            config,
            ignore_groups,
            ignore_unplaced,
            ignore_unlocalized
        )

    classify_ortho = \
        assign_clusters(
            classify_ortho,
            pprobs
        )

    if check_errors:
        col_index = pprobs.colindex
        row_index = pprobs.rowindex
        errors = ClusterErrorOrthogroups(samples=training_ortho.samples)
        for g in range(num(training_ortho.groups)):
            if classify_ortho.groups[g].cluster != training_ortho.groups[g].cluster:
                error = errors.new_group_copy(training_ortho.groups[g], append=True)
                for sample in training_ortho.samples:
                    error[sample.index] = training_ortho.groups[g][sample.index]
                error.cluster = (
                    training_ortho.groups[g].cluster,
                    classify_ortho.groups[g].cluster
                )
                error.probability = (
                    pprobs.matrix[row_index[training_ortho.groups[g].id]][col_index[training_ortho.groups[g].cluster]],
                    pprobs.matrix[row_index[classify_ortho.groups[g].id]][col_index[classify_ortho.groups[g].cluster]],
                )

        classify_ortho = errors

    return classify_ortho, pprobs
