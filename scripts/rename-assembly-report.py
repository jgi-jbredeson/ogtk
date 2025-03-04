#!/usr/bin/env python3

import os
import sys
import string

_ALPHA = str.maketrans({a:None for a in string.ascii_letters})
    
def format_number(field):
    number = field.translate(_ALPHA)
    length = len(number)
    if field == 'X' or field == 'Y' or \
       field == 'Z' or field == 'W':
        number = field
    elif length < 2:
        number = ('0' * (2-length)) + field
    else:
        number = field
    return number
        
def main(argv):
    assembly_report = open(argv[0], 'r')
    assembly_prefix = argv[1]

    unplaced_count = 0
    unlocalized_count = dict()
    for line in assembly_report:
        line   = line.strip()
        fields = line.split('\t')
        
        if line.startswith('#'):
            print(line)
            continue

        if fields[1] == 'assembled-molecule':
            if fields[3].lower().startswith('mito'):
                # Mitochondrion
                # Mitochondrial Plasmid
                fields[0] = assembly_prefix + fields[2]
            elif fields[3].lower().startswith('plastid') or \
                 fields[3].lower().startswith('chloro'):
                # Chloroplast
                # Plastid
                fields[0] = assembly_prefix + 'CP'
            else:  # Chromosome
                fields[0] = '%s%s' % (assembly_prefix, format_number(fields[2]))
                try:
                    unplaced_count = max(unplaced_count, int(fields[2]))
                except ValueError:
                    pass

        elif fields[1] == 'unlocalized-scaffold':
            if '_centromere_containing_' in fields[0]:
                fields[0] = '%sUn%s' % (assembly_prefix, format_number(fields[2]))
            else:
                fields[0] = '%sUn%s' % (assembly_prefix, format_number(fields[2]))

            if fields[0] in unlocalized_count:
                unlocalized_count[fields[0]] += 1
            else:
                unlocalized_count[fields[0]] = 1
            fields[0] += '.%d' % unlocalized_count[fields[0]]
                
        elif fields[1] == 'unplaced-scaffold':
            # # scaffold_
            # # Contig
            # # ChrUN-Ctg
            # # contig
            # # scaffold1A_scf_5
            # # super_
            # # unplaced
            # # Contig00058_ERROPOS234120
            # # Super_scaffold_552_2325899-2600059
            # # jcf7190000162957
            # # jcf7190000162957_0-116431
            # # jcf7190000000810_0-338178_0-118707
            # fields[0] = fields[0].lower()
            # if fields[0].startswith('super_scaffold_'):
            #     fields[0] = assembly_prefix + 'Un' + fields[0][15:]
            # elif fields[0].startswith('jcf7190000'):
            #     fields[0] = assembly_prefix + 'Un' + fields[0][10:]
            # elif fields[0].startswith('chrun-ctg'):
            #     fields[0] = assembly_prefix + 'Un' + fields[0][9:]
            # elif fields[0].startswith('scaffold_'):
            #     fields[0] = assembly_prefix + 'Un' + fields[0][9:]
            # elif fields[0].startswith('scaffold'):
            #     if '_scf_' in fields[0]:
            #         fields[0] = fields[0].replace('_scf_','')
            #     fields[0] = assembly_prefix + 'Un' + fields[0][8:]
            # elif fields[0].startswith('unplaced'):
            #     fields[0] = assembly_prefix + 'Un' + fields[0][8:]                
            # elif fields[0].startswith('super_'):
            #     fields[0] = assembly_prefix + 'Un' + fields[0][6:]
            # elif fields[0].startswith('contig'):
            #     if '_erropos' in fields[0]:
            #         fields[0] = fields[0].replace('_erropos','.')
            #     fields[0] = assembly_prefix + 'Un' + fields[0][6:]                
            # else:
            unplaced_count += 1
            fields[0] = '%sUn%d' % (assembly_prefix, unplaced_count)
                
                
        print('\t'.join(fields))


if __name__ == '__main__':
    main(sys.argv[1:])

    
