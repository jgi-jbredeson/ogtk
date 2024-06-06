#!/usr/bin/env python

import os
import sys
import re
import getopt

from math import inf as _POS_INF
from og.core.parsers.orthogroups import OrthoFinderOrthogroups
from og.constants import (
    _COLON,
    _COMMENT,
    _EMPTY,
    _EOL,
    _TAB,
)
__authors__ = 'Jessen V. Bredeson'
__program__ = os.path.basename(__file__)
__pkgname__ = '__PACKAGE_NAME__'
__version__ = '__PACKAGE_VERSION__'
__contact__ = '__PACKAGE_CONTACT__'
__purpose__ = 'Join OrthoFinder Orthogroups.tsv files'

num = len
max_missing_thresh = 2

_STRICT = 1


def format_intersection_counts_header():
    return "q_OG\tq_species\tq_members\tt_OG\tt_species\tt_members\tintersections\thits"



def format_intersection_counts_record(
        qry_id, qry_species, qry_members,
        trg_id, trg_species, trg_members, intersections, hits):
    return "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d" % (
        qry_id, qry_species, qry_members,
        trg_id, trg_species, trg_members,
        intersections, hits
    )



def get_setlist(orthogroup, loctable, unplaced_re=None, chronly=False):
    setlist = OrthoFinderOrthogroups()
    setlist.species = orthogroup.species
    setlist.ids = orthogroup.ids

    def _filter_un(setobj):
        return filter(lambda gname: not unplaced_re.match(loctable[species][gname][len(species):]), setobj)

    _filter = _filter_un if chronly and unplaced_re else lambda s: s
    
    for l in range(num(orthogroup.groups)):
        setlist.groups.append(set())
        for i, species in enumerate(orthogroup.species):
            if orthogroup.groups[l][i] is not None:
                setlist.groups[l].update(_filter(orthogroup.groups[l][i]))
    return setlist
    


def read_locus_bed(filename):
    records = dict()
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)
            
            assert num(fields) >= 6, "Too few BED fields, expected six: '%s'" % line

            fields[0] = fields[0].strip()
            fields[3] = fields[3].strip()
            
            if fields[3] in records:
                sys.stderr.write("Duplicate locus ID: %s" % fields[3])
            else:
                records[fields[3]] = fields[0]
                
    return records



def read_locus_bed_table(filename):
    locus_bed = dict()
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(maxsplit=1)
            fields[0] = fields[0].strip()
            fields[1] = fields[1].strip()
            locus_bed[fields[0]] = read_locus_bed(fields[1])
            
    return locus_bed



# def get_locus_species_map(ortho):
#     species_map = dict
#     for l in range(num(ortho.groups)):
#         for i in range(num(ortho.species)):
            

def _notNone(obj):
    return obj is not None and len(obj) > 0


def _map_locus_to_chr(locuslist, locusbed, species, unplaced_re, chronly=False):
    chrs = dict()
    for locusname in locuslist:
        if locusname not in locusbed:
            raise KeyError("Locus ID not in BED: %s" % locusname)

        chrname = locusbed[locusname]
        if chrname not in chrs:
            chrs[chrname] = 0        
        if chronly and \
           unplaced_re and \
           unplaced_re.match(chrname[len(species):]):
            continue
        chrs[chrname] += 1
            
    return chrs



def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else 'ERROR: %s\n\n' % message
    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage:   %s [options] <queryOG.tsv> <targetOG.tsv>\n" % __program__)
    stream.write("\n")
    stream.write("Options:\n")
    stream.write("  -a,--allow-multiple-best-targets\n")
    stream.write("     Allow queries to have multiple best (equal-scoring) targets. The\n")
    stream.write("     number of targets written per query is still restricted using `-Q`\n")
    stream.write("     By default, if the number of equally-scoring targets exceeds the\n")
    stream.write("     value of `-Q`, the query is left un-joined.\n")
    stream.write("\n")
    stream.write("  -b,--locus-bed-table <file>\n")
    stream.write("     Table of species ID and path to a BED file for each species. Used to\n")
    stream.write("     map locus IDs in the input orthogroups file to sequence names.\n")
    stream.write("\n")
    stream.write("  -c,--output-intersection-counts-file <file>\n")
    stream.write("     Write a counts table for intersetions between orthogroup sets.\n")
    stream.write("\n")
    stream.write("  -e,--regex-unplaced <regex>\n")
    stream.write("     Identify unplaced sequence using the specified regex.\n")
    stream.write("\n")
    stream.write("  -I,--ignore-unplaced-strictly\n")
    stream.write("     Strictly ignore unplaced sequences in filtering. If a cell in the\n")
    stream.write("     input orthogroups table contains no chromosomal sequences, that cell\n")
    stream.write("     then contains no members.\n")
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
    stream.write("  -n,--map-to-sequence-names\n")
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
    short_flags = 'hab:c:e:ILm:M:no:Q:s:S:T:'
    long_flags = (
        'help',
        'allow-multiple-best-targets',
        'locus-bed-table=',
        'regex-unplaced=',
        'ignore-unplaced-strictly',
        'left-outer-join',
        'min-intersecting-members=',
        'max-intersecting-members=',
        'min-intersecting-species=',
        'max-intersecting-species=',
        'max-intersecting-queries=',
        'max-intersecting-targets=',
        'map-to-sequence-names',
        'output-file=',
        'output-intersection-counts-file=',
    )
    try:
        options, arguments = getopt.getopt(argv, short_flags, long_flags)
    except getopt.GetoptError as error:
        usage(error)

    min_queries = 0
    max_queries = _POS_INF
    min_targets = 0
    max_targets = _POS_INF
    min_intersecting_members = 1
    max_intersecting_members = _POS_INF
    min_intersecting_species = 1
    max_intersecting_species = _POS_INF 
    locus_table = None
    left_outer_join = False
    write_seq_names = False
    output_file = sys.stdout
    intersection_counts_file = False
    ignore_unplaced = False
    regexp_unplaced = None
    allow_multiple_best_targets = False
    for flag, value in options:
        if   flag in ('-h','--help'):
            usage(exitcode=0)
        elif flag in ('-a','--allow-multiple-best-targets'):
            allow_multiple_best_targets = True
        elif flag in ('-e','--regex-unplaced'):
            regexp_unplaced = re.compile(value)
        elif flag in ('-I','--ignore-unplaced-strictly'):
            ignore_unplaced = _STRICT
        elif flag in ('-b','--locus-bed-table'):
            locus_table = value
        elif flag in ('-m','--min-intersecting-members'):
            min_intersecting_members = int(value)
        elif flag in ('-M','--max-intersecting-members'):
            max_intersecting_members = int(value)
        elif flag in ('-s','--min-intersecting-species'):
            min_intersecting_species = int(value)
        elif flag in ('-S','--max-intersecting-species'):
            max_intersecting_species = int(value)
        elif flag in ('-Q','--max-intersecting-queries'):
            max_queries = int(value) 
        elif flag in ('-T','--max-intersecting-targets'):
            max_targets = int(value)
        elif flag in ('-o','--output-file'):
            output_file = open(value, 'wt')
        elif flag in ('-n','--map-to-sequence-names'):
            write_seq_names = True
        elif flag in ('-c','--output-intersection-counts-file'):
            intersection_counts_file = open(value, 'wt')
        elif flag in ('-L','--left-outer-join'):
            left_outer_join = True
        
    if num(arguments) != 2:
        usage('Unexpected number of arguments')

    qry_ortho = OrthoFinderOrthogroups(arguments[0])
    trg_ortho = OrthoFinderOrthogroups(arguments[1])
    out_ortho = OrthoFinderOrthogroups()
    
    if locus_table is None:
        ignore_unplaced = False
        write_seq_names = False
    else:
        locus_table = read_locus_bed_table(locus_table)        
        for species in qry_ortho.species:
            if species not in locus_table:
                raise KeyError("Query species not found in locus BED files: %s" % species)
            
        for species in trg_ortho.species:
            if species not in locus_table:
                raise KeyError("Target species not found in locus BED files: %s" % species)

    qry_setlist = get_setlist(qry_ortho, locus_table, regexp_unplaced, ignore_unplaced)
    trg_setlist = get_setlist(trg_ortho, locus_table, regexp_unplaced, ignore_unplaced)

    qry_species = [ len(tuple(filter(_notNone, group))) for group in qry_ortho.groups ]
    trg_species = [ len(tuple(filter(_notNone, group))) for group in trg_ortho.groups ]
    add_species = sorted(set(trg_ortho.species) - set(qry_ortho.species))

    qry_species_index = dict(zip(qry_ortho.species, range(num(qry_ortho.species))))
    trg_species_index = dict(zip(trg_ortho.species, range(num(trg_ortho.species))))

    qry_counts = [0] * num(qry_ortho.groups)
    trg_counts = [0] * num(trg_ortho.groups)

    out_species = qry_ortho.species.copy()
    for species in trg_ortho.species:
        if species in add_species:
            out_species.append(species)

    if intersection_counts_file:
        intersection_counts_file.write(format_intersection_counts_header() + _EOL)

    output_file.write(out_ortho.format_orthogroups_header(out_species) + _EOL)

    intersection_list = []
    for q in range(num(qry_ortho.groups)):
        max_t = []
        max_num_intersecting_members = -1
        for t in range(num(trg_ortho.groups)):
            num_intersecting_members = num(qry_setlist.groups[q].intersection(trg_setlist.groups[t]))
            if num_intersecting_members > max_num_intersecting_members:
                max_num_intersecting_members = num_intersecting_members
                max_t = [t]
            elif num_intersecting_members == max_num_intersecting_members:
                max_t.append(t)

        intersection_list.extend(
            ((-max_num_intersecting_members, num(max_t), q, t) for t in reversed(max_t))
        )
                
    intersection_list.sort()

    
    for num_intersecting_members, num_t, q, t in intersection_list:
        if not allow_multiple_best_targets and num_t > max_queries:
            continue
        
        if (trg_counts[t] >= max_targets) or \
           (qry_counts[q] >= max_queries) or \
           (not (min_intersecting_members <= abs(num_intersecting_members) <= max_intersecting_members)):
            continue

        num_intersecting_species = 0
        intersecting_members = qry_setlist.groups[q].intersection(trg_setlist.groups[t])
        for species in qry_ortho.species:
            num_intersecting_species += num(
                set(qry_ortho.groups[q][qry_species_index[species]]).intersection(intersecting_members)
            )

        if not (min_intersecting_species <= num_intersecting_species <= max_intersecting_species):
            continue
            
        if intersection_counts_file:
            intersection_counts_file.write(
                format_intersection_counts_record(
                    qry_ortho.ids[q], qry_species[q], len(qry_setlist.groups[q]),
                    trg_ortho.ids[t], trg_species[t], len(trg_setlist.groups[t]),
                    abs(num_intersect), num_t
                ) + _EOL
            )
            
        out_id = qry_ortho.ids[q] + _COLON + trg_ortho.ids[t]
        out_group = qry_ortho.groups[q].copy()
        for species in add_species:
            out_group.append(trg_ortho.groups[t][trg_species_index[species]])

        if write_seq_names:
            for i in range(num(out_species)):
                out_group[i] = sorted(_map_locus_to_chr(
                    out_group[i],
                    locus_table[out_species[i]],
                    out_species[i],
                    regexp_unplaced,
                    chronly=ignore_unplaced
                ))
                            
        output_file.write(out_ortho.format_orthogroups_record(out_id, out_group) + _EOL)

        qry_counts[q] += 1
        trg_counts[t] += 1

    if left_outer_join:
        for q in range(num(qry_ortho.groups)):
            if qry_counts[q] > 0:
                continue
            if intersection_counts_file:
                intersection_counts_file.write(
                    format_intersection_counts_record(
                        qry_ortho.ids[q], qry_species[q], len(qry_setlist.groups[q]),
                        'NONE', 0, 0, 0, 0
                    ) + _EOL
                )
                
            out_id = qry_ortho.ids[q] + _COLON + 'NONE'
            out_group = qry_ortho.groups[q] + [tuple()] * num(add_species)
                
            if write_seq_names:
                for i in range(num(out_species)):
                    out_group[i] = sorted(_map_locus_to_chr(
                        out_group[i],
                        locus_table[out_species[i]],
                        out_species[i],
                        regexp_unplaced,
                        chronly=ignore_unplaced
                    ))
                    
            output_file.write(out_ortho.format_orthogroups_record(out_id, out_group) + _EOL)

    if intersection_counts_file:
        intersection_counts_file.close()
    output_file.close()
    
        
main(sys.argv[1:])
