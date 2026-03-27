
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
            "i_species\ti_members\t"
            "hits")


def format_intersection_counts_record(
        qry_id, qry_species, qry_members,
        trg_id, trg_species, trg_members, i_species, i_members, hits):
    return "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d" % (
        qry_id, qry_species, qry_members,
        trg_id, trg_species, trg_members,
        i_species, i_members, hits
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
                

def _isnotNone(obj):
    return obj is not None and len(obj) > 0


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
    stream.write("  -c,--output-intersection-counts-file <file>\n")
    stream.write("     Write a counts table for intersetions between orthogroup sets.\n")
    stream.write("\n")
    stream.write("  -I,--output-intersection-members\n")
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
    short_flags = 'haAc:ILm:M:No:Q:s:S:T:u'
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
        'output-intersection-members',
        'output-intersection-counts-file=',
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
    min_isec_members = 1
    max_isec_members = _POS_INF
    min_isec_species = 1
    max_isec_species = _POS_INF 
    left_outer_join = False
    output_seq_names = False
    output_all_members = False
    output_file = sys.stdout
    isec_counts_file = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-a','--allow-multiple-best-targets'):
            pass
        elif flag in ('-I','--output-all-members'):
            output_all_members = True
        elif flag in ('-m','--min-intersecting-members'):
            min_isec_members = int(value)
        elif flag in ('-M','--max-intersecting-members'):
            max_isec_members = int(value)
        elif flag in ('-s','--min-intersecting-species'):
            min_isec_species = int(value)
        elif flag in ('-S','--max-intersecting-species'):
            max_isec_species = int(value)
        elif flag in ('-Q','--max-intersecting-queries'):
            max_queries = int(value) 
        elif flag in ('-T','--max-intersecting-targets'):
            max_targets = int(value)
        elif flag in ('-o','--output-file'):
            output_file = open(value, 'wt')
        elif flag in ('-N','--output-sequence-names'):
            output_seq_names = True
        elif flag in ('-c','--output-intersection-counts-file'):
            isec_counts_file = open(value, 'wt')
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

    qry_species_count = [
        len(tuple(filter(_isnotNone, group))) for group in qry_ortho.groups
    ]
    trg_species_count = [
        len(tuple(filter(_isnotNone, group))) for group in trg_ortho.groups
    ]

    qry_members_index = index_members_by_ortho_ids(qry_ortho)
    trg_members_index = index_members_by_ortho_ids(trg_ortho)

    qry_species_index = index_list(qry_ortho.species)
    trg_species_index = index_list(trg_ortho.species)

    qry_ids_index = index_list(qry_ortho.ids)
    trg_ids_index = index_list(trg_ortho.ids)

    qry_mapped_counts = [0] * num(qry_ortho.groups)
    trg_mapped_counts = [0] * num(trg_ortho.groups)

    out_species = qry_ortho.species.copy()
    new_species = sorted(
        set(trg_ortho.species) - set(qry_ortho.species),
        key=trg_species_index.get
    )
    out_species.extend(new_species)

    
    if isec_counts_file:
        isec_counts_file.write(
            format_intersection_counts_header() + _EOL
        )

    output_file.write(
        out_ortho.format_orthogroups_header(out_species) + _EOL
    )

    isec_counts = dict()
    for member in qry_members_index:
        if member in trg_members_index:
            q = qry_members_index[member]
            t = trg_members_index[member]
            
            if (q,t) in isec_counts:
                isec_counts[(q,t)] += 1
            else:
                isec_counts[(q,t)] = 1

    qry_mapped_members = dict()
    trg_mapped_members = dict()
    for q,t in sorted(isec_counts, key=isec_counts.get, reverse=True):    
        num_isec_members = isec_counts[(q,t)]

        if (trg_mapped_counts[t] >= max_targets) or \
           (qry_mapped_counts[q] >= max_queries):
            continue

        num_isec_species = 0
        for s in range(num(qry_ortho.species)):
            num_isec_species += int(bool(
                set(qry_ortho.groups[q][s]) & trg_orthoset.groups[t]
            ))

        if not (min_isec_species <= num_isec_species <= max_isec_species):
            continue
        if not (min_isec_members <= num_isec_members <= max_isec_members):
            continue
        
        if isec_counts_file:
            isec_counts_file.write(
                format_intersection_counts_record(
                    qry_ortho.ids[q],
                    qry_species_count[q],
                    len(qry_orthoset.groups[q]),
                    trg_ortho.ids[t],
                    trg_species_count[t],
                    len(trg_orthoset.groups[t]),
                    num_isec_species,
                    num_isec_members,
                    0
                ) + _EOL
            )
            
        out_id = qry_ortho.ids[q] + _COLON + trg_ortho.ids[t]
        if output_all_members:
            out_group = [()] * num(qry_ortho.species)
            for s in range(num(qry_ortho.species)):
                out_group[s] = tuple(
                    set(qry_ortho.groups[q][s]) & trg_orthoset.groups[t]
                )
        else:
            out_group = qry_ortho.groups[q].copy()
            
        for species in new_species:
            out_group.append(
                trg_ortho.groups[t][trg_species_index[species]]
            )

        if left_outer_join and output_all_members:
            for member in (qry_orthoset.groups[q] & trg_orthoset.groups[t]):
                qry_members_index[member] = ~qry_members_index[member]
            
        if output_seq_names:
            for i in range(num(out_species)):
                out_group[i] = sorted(
                    map_loci_to_sequences(
                        out_group[i],
                        config.species[out_species[i]],
                        is_placed,
                        ignore_unplaced=False
                    )
                )
                            
        output_file.write(
            out_ortho.format_orthogroups_record(
                out_id, out_group
            ) + _EOL
        )

        qry_mapped_counts[q] += 1
        trg_mapped_counts[t] += 1

    if left_outer_join:
        for q in range(num(qry_ortho.groups)):
            if qry_mapped_counts[q] > 0:
                continue

            if output_all_members:
                num_members = 0
                for s in range(num(qry_ortho.species)):
                    if qry_ortho.groups[q][s] is None:
                        continue
                    for member in qry_ortho.groups[q][s]:
                        qry_members_index[member] = ~qry_members_index[member]
            
            if isec_counts_file:
                isec_counts_file.write(
                    format_intersection_counts_record(
                        qry_ortho.ids[q],
                        qry_species_count[q],
                        len(qry_orthoset.groups[q]),
                        'DISJOINT',
                        0,
                        0,
                        0,
                        0,
                        0
                    ) + _EOL
                )
                
            out_id = qry_ortho.ids[q] + _COLON + 'DISJOINT'
            out_group = qry_ortho.groups[q] + [()] * num(new_species)
                
            if output_seq_names:
                for s in range(num(out_species)):
                    out_group[s] = sorted(
                        map_loci_to_sequences(
                            out_group[s],
                            config.species[out_species[s]],
                            is_placed,
                            ignore_unplaced=False
                        )
                    )
                    
            output_file.write(
                out_ortho.format_orthogroups_record(
                    out_id, out_group
                ) + _EOL
            )
            
        if output_all_members:
            qq = set()
            for member in qry_members_index:
                if qry_members_index[member] < 0:
                    continue
                qq.add(qry_members_index[member])
                
            for q in sorted(qq):
                # if qry_mapped_counts[q] == 0:
                #     continue
                out_id = qry_ortho.ids[q] + _COLON + 'QRY_ONLY'
                out_group = [[] for s in range(num(out_species))]
                num_members = 0
                for s in range(num(qry_ortho.species)):
                    if qry_ortho.groups[q][s] is None:
                        continue
                    for member in qry_ortho.groups[q][s]:
                        if qry_members_index[member] < 0:
                            continue
                        out_group[s].append(member)
                        num_members += 1

                if num_members:
                    if isec_counts_file:
                        isec_counts_file.write(
                            format_intersection_counts_record(
                                qry_ortho.ids[q],
                                sum(map(bool, out_group)),
                                num_members,
                                'QRY_ONLY',
                                0,
                                0,
                                0,
                                0,
                                0
                            ) + _EOL
                        )
                    
                    if output_seq_names:
                        for i in range(num(out_species)):
                            out_group[i] = sorted(
                                map_loci_to_sequences(
                                    out_group[i],
                                    config.species[out_species[i]],
                                    is_placed,
                                    ignore_unplaced=False
                                )
                            )
                    output_file.write(
                        out_ortho.format_orthogroups_record(
                            out_id, out_group
                        ) + _EOL
                    )
            
    if isec_counts_file:
        isec_counts_file.close()
    output_file.close()
    
        
main(sys.argv[1:])
