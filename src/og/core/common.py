

_STRICT = 1
_LENIENT = 2

num = len

def map_loci_to_sequences(locus_names, namemap,
                          is_placed, ignore_unplaced=False):
    count = dict()
    for locus_name in locus_names:
        if locus_name not in namemap.loci:
            raise KeyError("Locus name not found: %s" % locus_name)

        increment_unit = 1
        sequence_name = namemap.loci[locus_name].chr
        if sequence_name not in namemap.references:
            raise KeyError("Sequence name not found: %s" % sequence_name)
        if is_placed(namemap.references[sequence_name]):
            sequence_name = \
                namemap.references[sequence_name].assigned_molecule
        elif ignore_unplaced:
            increment_unit = -1
            
        if sequence_name in count:
            count[sequence_name] += increment_unit
        else:
            count[sequence_name] = increment_unit
        
    return count


def filter_unplaced_sequences(sequence_names, namemap,
                              is_placed, ignore_unplaced=False):
    placed = set()
    unplaced = set()
    for sequence_name in sequence_names:
        if sequence_name not in namemap.references:
            raise KeyError("Sequence name not found: %s" % sequence_name)
        if ignore_unplaced and \
           not is_placed(namemap.references[sequence_name]):
            if ignore_unplaced == _LENIENT:
                unplaced.add(namemap.unplaced_id)
        else:
            placed.add(sequence_name)
    if ignore_unplaced == _LENIENT:
        if num(placed) < 1:
            placed = unplaced
    return placed

