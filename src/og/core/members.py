

_STRICT = 1
_LENIENT = 2

num = len

from og.core.parsers.assembly_report import is_placed as _placed


def int_placed(record, is_placed=_placed, ignore_unplaced=False):
    if ignore_unplaced and not is_placed(record):
        return -1 * int(ignore_unplaced == _LENIENT)
    return +1


def map_locus_to_sequence(locus_name, namemap, is_placed=_placed):
    if locus_name not in namemap.loci:
        raise KeyError(
            "Locus name not found in loci: " + str(locus_name)
        )
    sequence_name = namemap.loci[locus_name].chr
    if sequence_name not in namemap.references:
        raise KeyError(
            "Sequence name not found in references: " + str(sequence_name)
        )
    sequence_record = namemap.references[sequence_name]
    if is_placed(sequence_record):
        return sequence_record.assigned_molecule
    return sequence_name



def map_loci_to_sequences(locus_names, namemap,
                          is_placed=_placed, ignore_unplaced=False):
    sequences = dict()
    if not locus_names:
        return sequences

    for locus_name in locus_names:
        sequence_name = map_locus_to_sequence(locus_name, namemap, is_placed)
        increment_unit = int_placed(
            namemap.references[sequence_name],
            is_placed,
            ignore_unplaced
        )
        if increment_unit:
            if increment_unit < 0:
                sequence_name = namemap.unplaced_id
        else:
            sequence_name = None
        sequences[locus_name] = sequence_name

    return sequences
        

        
def map_loci_to_sequence_counts(locus_names, namemap,
                                is_placed=_placed, ignore_unplaced=False):
    count = dict()
    if not locus_names:
        return count
    
    for locus_name in locus_names:
        sequence_name = map_locus_to_sequence(locus_name, namemap, is_placed)
        increment_unit = int_placed(
            namemap.references[sequence_name],
            is_placed,
            ignore_unplaced
        )

        if increment_unit:
            try:
                count[sequence_name] += increment_unit
            except KeyError:
                count[sequence_name] = increment_unit
        
    return count


def filter_unplaced_sequences(sequence_names, namemap, is_placed=_placed,
                              ignore_unplaced=False, aggregate_unplaced=True):
    placed = set()
    unplaced = set()
    if not sequence_names:
        return placed
    
    for sequence_name in sequence_names:
        if sequence_name not in namemap.references:
            raise KeyError(
                "Sequence name not found in references: %s" % sequence_name
            )

        increment_unit = 1
        if is_placed(namemap.references[sequence_name]):
            sequence_name = \
                namemap.references[sequence_name].assigned_molecule
        elif ignore_unplaced:
            if ignore_unplaced == _LENIENT:
                unplaced.add(
                    namemap.unplaced_id \
                    if   aggregate_unplaced \
                    else sequence_name
                )
            increment_unit = 0

        if increment_unit:
            placed.add(sequence_name)
            
    if ignore_unplaced == _LENIENT:
        if num(placed) < 1:
            placed = unplaced

    return placed

