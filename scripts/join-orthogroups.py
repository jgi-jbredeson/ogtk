
import os
import sys
import getopt

from math import inf as _POS_INF
from og.core.utils import index_list
from og.core.members import map_loci_to_sequence_counts
from og.core.parsers.config import SampleConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
from og.constants import (
    _COMMENT,
    _EMPTY,
    _EOL
)

__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Join OrthoFinder Orthogroups.tsv files'

num = len
max_missing_thresh = 2


def format_id(qry_id, trg_id, sep=':', reverse=False):
    if reverse:
        return str(trg_id) + sep + str(qry_id)
    else:
        return str(qry_id) + sep + str(trg_id)


def format_intersection_counts_header():
    return ("q_OG\tq_samples\tq_members\t"
            "t_OG\tt_samples\tt_members\t"
            "i_samples\ti_members")


def format_intersection_counts_record(
        qry_id, qry_samples, qry_members,
        trg_id, trg_samples, trg_members, i_samples, i_members):
    return "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d" % (
        qry_id, qry_samples, qry_members,
        trg_id, trg_samples, trg_members,
        i_samples, i_members
    )


def index_members_by_ortho_ids(ortho):
    index = dict()
    for g in range(num(ortho.groups)):
        for s in range(num(ortho.samples)):
            if not ortho.groups[g][s]:
                continue
            for member in ortho.groups[g][s]:
                index[member] = g
    return index


def get_orthoset(ortho):
    orthoset = OrthoFinderOrthogroups(samples=["set"])
    for g in range(num(ortho.groups)):
        group = orthoset.new_group(append=True)
        group.id = ortho.groups[g].id
        for s in range(num(ortho.samples)):
            if not ortho.groups[g][s]:
                continue
            group[0].update(ortho.groups[g][s])
                
    return orthoset
                


def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <queryOG.tsv> <targetOG.tsv> [in.yaml]\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    # stream.write("  -a,--allow-multiple-best-targets\n")
    # stream.write("     For each query, write equivalently-scoring targets.\n")
    # stream.write("     This option modifies `--max-intersecting-queries`\n")
    # stream.write("     behavior.\n")
    # stream.write("\n")
    stream.write("  -c,--output-intersecting-counts-file <file>\n")
    stream.write("     Write a counts table for intersetions between orthogroup sets.\n")
    stream.write("\n")
    stream.write("  -I,--output-intersecting-members\n")
    stream.write("     For each query-target pair, output only intersecting members [OG]\n")
    stream.write("\n")
    stream.write("  -L,--left-outer-join\n")
    stream.write("     Perform a left outer join [inner]\n")
    stream.write("\n")
    stream.write("  -m,--min-intersecting-members <uint>\n")
    stream.write("     Minimum number of intersecting members permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -M,--max-intersecting-members <uint>\n")
    stream.write("     Maximum number of intersecting members permitted per orthogroup [inf]\n")
    stream.write("\n")
    stream.write("  -N,--output-sequence-names\n")
    stream.write("     Write sequence names to the output orthogroups table [locus IDs]\n")
    stream.write("\n")
    stream.write("  -o,--output-file <file>\n")
    stream.write("     Write output to file [stdout]\n")
    stream.write("\n")
    stream.write("  -r,--reverse-id-order\n")
    stream.write("     By default, the IDs of joined orthogroups are concatenated together\n")
    stream.write("     as `qryOGID:trgOGID`, this option instead outputs `trgOGID:qryOGID`\n")
    stream.write("\n")
    stream.write("  -s,--min-intersecting-samples <uint>\n")
    stream.write("     Minimum number of intersecting samples permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -S,--max-intersecting-samples <uint>\n")
    stream.write("     Maximum number of intersecting samples permitted per orthogroup [inf]\n")
    stream.write("\n")
    stream.write("  -Q,--max-intersecting-queries <uint>\n")
    stream.write("     Maximum number of equally-scoring target orthogroups per query [inf]\n")
    stream.write("\n")
    stream.write("  -T,--max-intersecting-targets <uint>\n")
    stream.write("     Maximum number of intersecting target orthogroups counted across all\n")
    stream.write("     queries [inf]\n")
    stream.write("\n")
    stream.write("  -u,--ignore-unlocalized\n")
    stream.write("     Map the names of loci on (placed but) unlocalized sequences to their\n")
    stream.write("     designated sequence names, not to their placed chromosome names.\n")
    stream.write("     (only effective when --output-sequence-names is enabled)\n")
    stream.write("\n")
    stream.write("  -x,--output-intersecting-strict\n")
    stream.write("     When performing the intersection, drop individual members from the\n")
    stream.write("     query that are not also in the target\n")
    stream.write("\n")
    stream.write("  -h,--help\n")
    stream.write("     Print this help message and exit\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("\n")
    stream.write("Notes:\n")
    stream.write("  1. queryOG.tsv and targetOG.tsv files are Orthogroups.tsv files with\n")
    stream.write("     headers defined.\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    sys.exit(exitcode)

    
def main(argv):
    short_flags = 'haAc:ILm:M:No:Q:rs:S:T:ux'
    long_flags = (
        'help',
        'left-outer-join',
        'min-intersecting-members=',
        'max-intersecting-members=',
        'min-intersecting-samples=','min-intersecting-species=',
        'max-intersecting-samples=','max-intersecting-species=',
        'max-intersecting-queries=',
        'max-intersecting-targets=',
        'output-sequence-names',
        'output-file=',
        'allow-multiple-best-targets',
        'output-intersecting-members',
        'output-intersecting-strict',
        'output-intersecting-counts-file=',
        'ignore-unlocalized',
        'reverse-id-order',
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    is_placed = _placed
    min_queries = 0
    max_queries = _POS_INF
    min_targets = 0
    max_targets = _POS_INF
    reverse_id_order = False
    min_intersect_members = 1
    max_intersect_members = _POS_INF
    min_intersect_samples = 1
    max_intersect_samples = _POS_INF
    left_outer_join = False
    output_seq_names = False
    output_all_best = False
    output_all_members = False
    output_file = sys.stdout
    intersect_counts_file = False
    intersect_strict = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-I','--output-all-members'):
            output_all_members = True
        elif flag in ('-m','--min-intersecting-members'):
            min_intersect_members = int(value)
        elif flag in ('-M','--max-intersecting-members'):
            max_intersect_members = int(value)
        elif flag in ('-s','--min-intersecting-samples','--min-intersecting-species'):
            min_intersect_samples = int(value)
        elif flag in ('-S','--max-intersecting-samples','--max-intersecting-species'):
            max_intersect_samples = int(value)
        elif flag in ('-r','--reverse-id-order'):
            reverse_id_order = True
        elif flag in ('-Q','--max-intersecting-queries'):
            max_queries = int(value) 
        elif flag in ('-T','--max-intersecting-targets'):
            max_targets = int(value)
        elif flag in ('-o','--output-file'):
            output_file = open(value, 'wt')
        elif flag in ('-N','--output-sequence-names'):
            output_seq_names = True
        elif flag in ('-c','--output-intersecting-counts-file'):
            intersect_counts_file = open(value, 'wt')
        elif flag in ('-x','--output-intersecting-strict'):
            intersect_strict = True
        elif flag in ('-L','--left-outer-join'):
            left_outer_join = True
        elif flag in ('-u','--ignore-unlocalized'):
            is_placed = _localized

    if num(arguments) < 2 or num(arguments) > 3:
        usage('Unexpected number of arguments')

    qry_ortho = OrthoFinderOrthogroups(arguments[0])
    trg_ortho = OrthoFinderOrthogroups(arguments[1])
        
    if output_seq_names:
        if num(arguments) != 3:
            usage('`--output-sequence-names` detected, but no YAML file given')
        config = SampleConfigFile(arguments[2], load_files=True, map_assigned_molecule=True)

        for sample in qry_ortho.samples:
            if sample.id not in config.samples:
                raise KeyError(
                    "Sample not found in YAML file: '%s'" % str(sample.id)
                )
        for sample in trg_ortho.samples:
            if sample.id not in config.samples:
                raise KeyError(
                    "Sample not found in YAML file: '%s'" % str(sample.id)
                )
    else:
        config = None
    
    qry_orthoset = get_orthoset(qry_ortho)
    trg_orthoset = get_orthoset(trg_ortho)

    qry_member_indices = index_members_by_ortho_ids(qry_ortho)
    trg_member_indices = index_members_by_ortho_ids(trg_ortho)

    qry_sample_indices = {s.id:s.index for s in qry_ortho.samples}
    trg_sample_indices = {s.id:s.index for s in trg_ortho.samples}
    
    qry_mapped_counts = [0] * num(qry_ortho.groups)
    trg_mapped_counts = [0] * num(trg_ortho.groups)

    qry_samples = [s.id for s in qry_ortho.samples]
    trg_samples = sorted(
        set(s.id for s in trg_ortho.samples) - set(qry_samples),
        key=trg_sample_indices.get
    )
    mrg_ortho = OrthoFinderOrthogroups(samples=(qry_samples + trg_samples))
    mrg_sample_indices = {s.id:s.index for s in mrg_ortho.samples}

    intersect_counts = dict()
    for member in qry_member_indices:
        if member in trg_member_indices:
            q = qry_member_indices[member]
            t = trg_member_indices[member]
            
            if (q,t) in intersect_counts:
                intersect_counts[(q,t)] += 1
            else:
                intersect_counts[(q,t)] = 1

    mrg_index = 0
    out_counts = list()
    qry_only_member_dict = dict()
    for q,t in sorted(intersect_counts, key=intersect_counts.get, reverse=True):
        if ((trg_mapped_counts[t] >= max_targets) or
            (qry_mapped_counts[q] >= max_queries)):
            continue
        
        num_intersect_members = intersect_counts[(q,t)]
        num_intersect_samples = sum(
            bool(qry_ortho.groups[q][s] & trg_orthoset.groups[t][0]) \
            for s in range(num(qry_samples))
        )

        
        if not (min_intersect_samples <= num_intersect_samples <= max_intersect_samples):
            continue
        if not (min_intersect_members <= num_intersect_members <= max_intersect_members):
            continue

        mrg_group = mrg_ortho.new_group(append=True)
        mrg_group.id = format_id(
            qry_ortho.groups[q].id,
            trg_ortho.groups[t].id,
            reverse=reverse_id_order
        )
        if output_all_members:
            for sample_id in qry_samples:
                qs = qry_sample_indices[sample_id]
                ms = mrg_sample_indices[sample_id]
                mrg_group[ms] = qry_ortho.groups[q][qs] & trg_orthoset.groups[t][0]
                
                if not intersect_strict:
                    for member in (qry_ortho.groups[q][qs] - mrg_group[ms]):
                        if member in qry_only_member_dict:
                            qry_only_member_dict[member].append((mrg_index, ms))
                        else:
                            qry_only_member_dict[member] = [(mrg_index, ms)]
                        
        else:
            for sample_id in qry_samples:
                qs = qry_sample_indices[sample_id]
                ms = mrg_sample_indices[sample_id]
                mrg_group[ms] = qry_ortho.groups[q][qs]
                
        for sample_id in trg_samples:
            ts = trg_sample_indices[sample_id]
            ms = mrg_sample_indices[sample_id]
            mrg_group[ms] = trg_ortho.groups[t][ts]

        for sample_id in qry_samples:
            ms = mrg_sample_indices[sample_id]
            for member in mrg_group[ms]:
                if qry_member_indices[member] >= 0:
                    qry_member_indices[member] = \
                        ~qry_member_indices[member]
            
        out_counts.append((
            qry_ortho.groups[q].id,
            qry_ortho.groups[q].num_samples,
            qry_ortho.groups[q].num_members,
            trg_ortho.groups[t].id,
            trg_ortho.groups[t].num_samples,
            trg_ortho.groups[t].num_members,
            num_intersect_samples,
            num_intersect_members,
        ))
        
        qry_mapped_counts[q] += 1
        trg_mapped_counts[t] += 1
        mrg_index += 1

    # if dropped members are uniquely mapped, add them back in by back-filling:
    for member in qry_only_member_dict:
        if qry_member_indices[member] >= 0 and \
           num(qry_only_member_dict[member]) == 1:
            g,s = qry_only_member_dict[member].pop()
            mrg_ortho.groups[g][s].add(member)
            qry_member_indices[member] = \
                ~qry_member_indices[member]
    
    qry_only_member_dict = None

    if left_outer_join:
        for q in range(num(qry_ortho.groups)):
            if qry_mapped_counts[q] > 0:
                continue

            if output_all_members:
                for s in range(num(qry_samples)):
                    for member in qry_ortho.groups[q][s]:
                        if qry_member_indices[member] >= 0:
                            qry_member_indices[member] = \
                                ~qry_member_indices[member]

            mrg_group = mrg_ortho.new_group(append=True)
            mrg_group.id = format_id(
                qry_ortho.groups[q].id,
                'DISJOINT',
                reverse=reverse_id_order
            )
            for sample_id in qry_samples:
                qs = qry_sample_indices[sample_id]
                ms = mrg_sample_indices[sample_id]
                mrg_group[ms] = qry_ortho.groups[q][qs]
            
            out_counts.append((
                qry_ortho.groups[q].id,
                qry_ortho.groups[q].num_samples,
                qry_ortho.groups[q].num_members,
                'DISJOINT',
                0,
                0,
                qry_ortho.groups[q].num_samples,
                qry_ortho.groups[q].num_members,
            ))
            mrg_index += 1
                    
        if output_all_members:
            for q in set(filter((0).__le__, qry_member_indices.values())):
                mrg_group = mrg_ortho.new_group(append=False)
                mrg_group.id = format_id(
                    qry_ortho.groups[q].id,
                    'QRY_ONLY',
                    reverse=reverse_id_order
                )
                for sample_id in qry_samples:
                    qs = qry_sample_indices[sample_id]
                    ms = mrg_sample_indices[sample_id]
                    for member in qry_ortho.groups[q][qs]:
                        if qry_member_indices[member] >= 0:
                            mrg_group[ms].add(member)

                if mrg_group.num_members:
                    mrg_ortho.groups.append(mrg_group)
                    out_counts.append((
                        qry_ortho.groups[q].id,
                        qry_ortho.groups[q].num_samples,
                        qry_ortho.groups[q].num_members,
                        'QRY_ONLY',
                        0,
                        0,
                        mrg_group.num_samples,
                        mrg_group.num_members,
                    ))

    if output_seq_names:
        for group in mrg_ortho.groups:
            for sample in mrg_ortho.samples:
                group[sample.index] = \
                    sorted(map_loci_to_sequence_counts(
                        group[sample.index],
                        config.samples[sample.id],
                        is_placed,
                        ignore_unplaced=False
                    )
                )

    mrg_ortho.to_file(output_file)
    
    if intersect_counts_file:
        intersect_counts_file.write(
            format_intersection_counts_header() + _EOL
        )
        for record in out_counts:
            intersect_counts_file.write(
                format_intersection_counts_record(*record) + _EOL
            )
        intersect_counts_file.close()



main(sys.argv[1:])
