

def map_loci_to_sequences(locus_names, namemap, is_placed, ignore_unplaced=False):
    count = dict()
    for locus_name in locus_names:
        if locus_name not in namemap.loci:
            raise KeyError("Locus name not found: %s" % locus_name)

        increment_unit = 1
        sequence_name = namemap.loci[locus_name].chr
        if is_placed(namemap.references[sequence_name]):
            sequence_name = namemap.references[sequence_name].assigned_molecule

        elif ignore_unplaced:
            increment_unit = -1
            
        if sequence_name in count:
            count[sequence_name] += increment_unit
        else:
            count[sequence_name] = increment_unit
        
    return count
