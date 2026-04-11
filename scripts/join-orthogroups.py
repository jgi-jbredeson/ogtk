
import os
import sys
import getopt

from math import inf as _POS_INF
from og.core.utils import index_list
from og.core.members import map_loci_to_sequences
from og.core.parsers.config import SpeciesConfigFile
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.core.parsers.assembly_report import is_chr as _localized
from og.core.parsers.assembly_report import is_placed as _placed
from og.constants import (
    _COLON,
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


def format_intersection_counts_header():
    return ("q_OG\tq_species\tq_members\t"
            "t_OG\tt_species\tt_members\t"
            "i_species\ti_members")


def format_intersection_counts_record(
        qry_id, qry_species, qry_members,
        trg_id, trg_species, trg_members, i_species, i_members):
    return "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d" % (
        qry_id, qry_species, qry_members,
        trg_id, trg_species, trg_members,
        i_species, i_members
    )


def index_members_by_ortho_ids(ortho):
    index = dict()
    for g in range(num(ortho.groups)):
        for s in range(num(ortho.species)):
            if ortho.groups[g][s] is None:
                continue
            for member in ortho.groups[g][s]:
                index[member] = g
    return index


def get_orthoset(orthogroup):
    orthoset = OrthoFinderOrthogroups()
    orthoset.species = orthogroup.species
    orthoset.ids = orthogroup.ids

    for l in range(num(orthogroup.groups)):
        orthoset.groups.append(set())
        for i, species in enumerate(orthogroup.species):
            if orthogroup.groups[l][i] is not None:
                orthoset.groups[l].update(orthogroup.groups[l][i])
                
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
    stream.write("  -s,--min-intersecting-species <uint>\n")
    stream.write("     Minimum number of intersecting species permitted per orthogroup [1]\n")
    stream.write("\n")
    stream.write("  -S,--max-intersecting-species <uint>\n")
    stream.write("     Maximum number of intersecting species permitted per orthogroup [inf]\n")
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
    short_flags = 'haAc:ILm:M:No:Q:s:S:T:ux'
    long_flags = (
        'help',
        'left-outer-join',
        'min-intersecting-members=',
        'max-intersecting-members=',
        'min-intersecting-species=',
        'max-intersecting-species=',
        'max-intersecting-queries=',
        'max-intersecting-targets=',
        'output-sequence-names',
        'output-file=',
        'allow-multiple-best-targets',
        'output-intersecting-members',
        'output-intersecting-strict',
        'output-intersecting-counts-file=',
        'ignore-unlocalized'
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
    min_intersect_members = 1
    max_intersect_members = _POS_INF
    min_intersect_species = 1
    max_intersect_species = _POS_INF 
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
        elif flag in ('-s','--min-intersecting-species'):
            min_intersect_species = int(value)
        elif flag in ('-S','--max-intersecting-species'):
            max_intersect_species = int(value)
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
    out_ortho = OrthoFinderOrthogroups()
        
    if output_seq_names:
        if num(arguments) != 3:
            usage('--output-sequence-names requested, but no YAML file given')
        config = SpeciesConfigFile(arguments[2], load_files=True, map_assigned_molecule=True)

        for species_id in qry_ortho.species:
            if species_id not in config.species:
                raise KeyError("Species not found in YAML file: '%s'" % (
                    str(species_id)
                ))
        for species_id in trg_ortho.species:
            if species_id not in config.species:
                raise KeyError("Species not found in YAML file: '%s'" % (
                    str(species_id)
                ))
    else:
        config = None
    
    qry_orthoset = get_orthoset(qry_ortho)
    trg_orthoset = get_orthoset(trg_ortho)

    qry_members_indices = index_members_by_ortho_ids(qry_ortho)
    trg_members_indices = index_members_by_ortho_ids(trg_ortho)
    
    qry_species_indices = index_list(qry_ortho.species)
    trg_species_indices = index_list(trg_ortho.species)
    qry_species_counts = [sum(map(bool, g)) for g in qry_ortho.groups]
    trg_species_counts = [sum(map(bool, g)) for g in trg_ortho.groups]
    
    qry_ids_indices = index_list(qry_ortho.ids)
    trg_ids_indices = index_list(trg_ortho.ids)

    qry_mapped_members = [0] * num(qry_ortho.groups)
    trg_mapped_members = [0] * num(trg_ortho.groups)    
    qry_mapped_counts = [0] * num(qry_ortho.groups)
    trg_mapped_counts = [0] * num(trg_ortho.groups)

    
    out_ortho.species = qry_ortho.species.copy()
    new_species = sorted(
        set(trg_ortho.species) - set(qry_ortho.species),
        key=trg_species_indices.get
    )
    out_ortho.species.extend(new_species)
    
    intersect_counts = dict()
    for member in qry_members_indices:
        if member in trg_members_indices:
            q = qry_members_indices[member]
            t = trg_members_indices[member]
            
            if (q,t) in intersect_counts:
                intersect_counts[(q,t)] += 1
            else:
                intersect_counts[(q,t)] = 1

    out_index = 0
    out_counts = list()
    qry_only_members_dict = dict()
    for q,t in sorted(intersect_counts, key=intersect_counts.get, reverse=True):
        if (trg_mapped_counts[t] >= max_targets) or \
           (qry_mapped_counts[q] >= max_queries):
            continue
        
        num_intersect_members = intersect_counts[(q,t)]
        num_intersect_species = sum([
            int(bool(set(qry_ortho.groups[q][s]) & trg_orthoset.groups[t])) \
            for s in range(num(qry_ortho.species))
        ])

        if not (min_intersect_species <= num_intersect_species <= max_intersect_species):
            continue
        if not (min_intersect_members <= num_intersect_members <= max_intersect_members):
            continue

        if output_all_members:
            out_group = [[] for s in range(num(qry_ortho.species))]
            for s in range(num(qry_ortho.species)):
                qry_set = set(qry_ortho.groups[q][s])
                
                out_group[s] = list(qry_set & trg_orthoset.groups[t])
                if not intersect_strict:
                    for member in (qry_set - trg_orthoset.groups[t]):
                        if member not in qry_only_members_dict:
                            qry_only_members_dict[member] = []
                        qry_only_members_dict[member].append((out_index,s))
        else:
            out_group = qry_ortho.groups[q].copy()
                
        for species in new_species:
            out_group.append(
                trg_ortho.groups[t][trg_species_indices[species]]
            )

        for s in range(num(qry_ortho.species)):
            if out_group[s]:
                for member in out_group[s]:
                    if qry_members_indices[member] >= 0:
                        qry_members_indices[member] = \
                            ~qry_members_indices[member]
            
        out_ortho.groups.append(
            out_group
        )
        out_ortho.ids.append(
            qry_ortho.ids[q] + _COLON + trg_ortho.ids[t]
        )
        out_counts.append((
            qry_ortho.ids[q],
            qry_species_counts[q],
            num(qry_orthoset.groups[q]),
            trg_ortho.ids[t],
            trg_species_counts[t],
            num(trg_orthoset.groups[t]),
            num_intersect_species,
            num_intersect_members,
        ))
        
        qry_mapped_members[q] = num_intersect_members
        qry_mapped_counts[q] += 1
        trg_mapped_counts[t] += 1
        out_index += 1

    # if dropped members are uniquely mapped, add them back in by back-filling:
    for member in qry_only_members_dict:
        if qry_members_indices[member] >= 0 and \
           num(qry_only_members_dict[member]) == 1:
            g,s = qry_only_members_dict[member].pop()
            out_ortho.groups[g][s].append(member)
            qry_members_indices[member] = \
                ~qry_members_indices[member]
    
    qry_only_members_dict = None

    if left_outer_join:
        for q in range(num(qry_ortho.groups)):
            if qry_mapped_counts[q] > 0:
                continue

            if output_all_members:
                num_members = 0
                for s in range(num(qry_ortho.species)):
                    if qry_ortho.groups[q][s]:
                        for member in qry_ortho.groups[q][s]:
                            if qry_members_indices[member] >= 0:
                                qry_members_indices[member] = \
                                    ~qry_members_indices[member]
                
            out_ortho.groups.append(
                qry_ortho.groups[q] + [()] * num(new_species)
            )
            out_ortho.ids.append(
                qry_ortho.ids[q] + _COLON + 'DISJOINT'
            )
            out_counts.append((
                qry_ortho.ids[q],
                qry_species_counts[q],
                num(qry_orthoset.groups[q]),
                'DISJOINT',
                0,
                0,
                qry_species_counts[q],
                num(qry_orthoset.groups[q]),                
            ))

                    
        if output_all_members:
            for q in set(filter((0).__le__, qry_members_indices.values())):
                num_members = 0
                out_group = [[] for s in range(num(out_ortho.species))]
                for s in range(num(qry_ortho.species)):
                    if qry_ortho.groups[q][s]:
                        for member in qry_ortho.groups[q][s]:
                            if qry_members_indices[member] >= 0:
                                out_group[s].append(member)
                                num_members += 1
                if num_members:
                    out_ortho.groups.append(
                        out_group
                    )
                    out_ortho.ids.append(
                        qry_ortho.ids[q] + _COLON + 'QRY_ONLY'
                    )
                    out_counts.append((
                        qry_ortho.ids[q],
                        qry_species_counts[q],
                        num(qry_orthoset.groups[q]),
                        'QRY_ONLY',
                        0,
                        0,
                        sum(map(bool, out_group)),
                        num_members,                        
                    ))

    if output_seq_names:
        for g in range(num(out_ortho.groups)):
            for s in range(num(out_ortho.species)):
                out_ortho.groups[g][s] = sorted(
                    map_loci_to_sequences(
                        out_ortho.group[g][s],
                        config.species[out_ortho.species[s]],
                        is_placed,
                        ignore_unplaced=False
                    )
                )

    out_ortho.to_table(file=output_file)
    
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
