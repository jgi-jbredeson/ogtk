
_LENIENT = 1
_STRICT = 2


num = len

from og.core.assembly_report import is_placed as _placed


def int_placed(record, is_placed=_placed, ignore_unplaced=False):
    if ignore_unplaced and not is_placed(record):
        return -1 * int(ignore_unplaced == _LENIENT)
    return +1


def map_locus_to_reference(locus_name, namemap, is_placed=False):
    if locus_name not in namemap.loci:
        raise KeyError(
            "Locus name not found in loci: " + str(locus_name)
        )
    if not is_placed:
        is_placed = _placed
    reference_name = namemap.loci[locus_name].chr
    if reference_name not in namemap.references:
        raise KeyError(
            "Reference name not found in references: " + str(reference_name)
        )
    reference_record = namemap.references[reference_name]
    if is_placed(reference_record):
        return reference_record.assigned_molecule
    return reference_name



def map_loci_to_references(locus_names, namemap,
                           is_placed=False, ignore_unplaced=False):
    references = dict()
    if not locus_names:
        return references
    if not is_placed:
        is_placed = _placed
    for locus_name in locus_names:
        reference_name = map_locus_to_reference(locus_name, namemap, is_placed)
        increment_unit = int_placed(
            namemap.references[reference_name],
            is_placed,
            ignore_unplaced
        )
        if increment_unit:
            if increment_unit < 0:
                reference_name = namemap.unplaced_id
        else:
            reference_name = None
        references[locus_name] = reference_name

    return references
        

        
def map_loci_to_reference_counts(locus_names, namemap,
                                is_placed=False, ignore_unplaced=False):
    count = dict()
    if not locus_names:
        return count
    if not is_placed:
        is_placed = _placed
    for locus_name in locus_names:
        reference_name = map_locus_to_reference(locus_name, namemap, is_placed)
        increment_unit = int_placed(
            namemap.references[reference_name],
            is_placed,
            ignore_unplaced
        )

        if increment_unit:
            try:
                count[reference_name] += increment_unit
            except KeyError:
                count[reference_name] = increment_unit
        
    return count


def filter_unplaced_references(reference_names, namemap, is_placed=False,
                               ignore_unplaced=False, aggregate_unplaced=True):
    placed = set()
    unplaced = set()
    if not reference_names:
        return placed
    if not is_placed:
        is_placed = _placed
    for reference_name in reference_names:
        if reference_name not in namemap.references:
            raise KeyError(
                "Reference name not found in references: %s" % reference_name
            )

        increment_unit = 1
        if is_placed(namemap.references[reference_name]):
            reference_name = \
                namemap.references[reference_name].assigned_molecule
        elif ignore_unplaced:
            if ignore_unplaced == _LENIENT:
                unplaced.add(
                    namemap.unplaced_id \
                    if   aggregate_unplaced \
                    else reference_name
                )
            increment_unit = 0

        if increment_unit:
            placed.add(reference_name)
            
    if ignore_unplaced == _LENIENT:
        if num(placed) < 1:
            placed = unplaced

    return placed


def map_members_to_group_index(ortho, samples=()):
    index = dict()
    samples = samples or ortho.samples
    for g, group in enumerate(ortho.groups):
        for sample in samples:
            if group[sample.index]:
                for member in group[sample.index]:
                    index[member] = g
    return index
