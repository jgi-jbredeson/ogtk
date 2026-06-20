
from og.core.generic import Enumerative

ENUM_OUTPUT = Enumerative(
    ANCHORS='A',
    ORTHOFINDER='F',
    ORTHOVENN='V'
)
ENUM_SORTBY = Enumerative(
    UNSORTED=0x0,
    OGID=0x1,
    MEMBERS=0x2,
    SAMPLES=0x4
)
ENUM_CALLKEY = Enumerative(
    OUTPUT_REFS='output-reference-names',
    INPUT_REFS='input-reference-names',
    PREFIX_SAMPLES='prefix-samples',
    REMOVE_SAMPLES='remove-samples',
    SUBSET_SAMPLES='subset-samples',
    SORTBY='sort-by'
)

num = len



def view(
        ortho,
        config=None,
        call_order=(),
        prefix_delim='|',
        replace_delim='|',
        ignore_unplaced=None,
        ignore_unlocalized=False,
        input_reference_names=False,
        output_reference_names=False,
):
    call_function = dict()
    def call_output_reference_names(*args):
        ortho.map_loci_to_references(
            config,
            is_placed=ignore_unlocalized,
            ignore_unplaced=ignore_unplaced,
            inplace=True
        )
    call_function[ENUM_CALLKEY.OUTPUT_REFS] = call_output_reference_names
            
    def call_input_reference_names(*args):
        ortho.filter_unplaced_references(
            config,
            is_placed=ignore_unlocalized,
            ignore_unplaced=ignore_unplaced,
            aggregate_unplaced=True,
            inplace=True
        )
    call_function[ENUM_CALLKEY.INPUT_REFS] = call_input_reference_names

    def call_prefix_samples(*args):
        ortho.prepend_sample_prefixes(
            prefix_delim,
            replace_delim,
            inplace=True
        )
    call_function[ENUM_CALLKEY.PREFIX_SAMPLES] = call_prefix_samples
    
    def call_remove_samples(*args):
        ortho.remove_sample_prefixes(
            prefix_delim,
            inplace=True
        )
    call_function[ENUM_CALLKEY.REMOVE_SAMPLES] = call_remove_samples
    
    def call_sortby(sortby):
        if abs(sortby) & ENUM_SORTBY.OGID:
            key = lambda g: g.id
        elif abs(sortby) == (ENUM_SORTBY.SAMPLES | ENUM_SORTBY.MEMBERS):
            key = lambda g: (g.num_samples, g.num_members)
        elif abs(sortby) & ENUM_SORTBY.MEMBERS:
            key = lambda g: g.num_members
        elif abs(sortby) & ENUM_SORTBY.SAMPLES:
            key = lambda g: g.num_samples
        else:
            raise AssertionError('Invalid sort value: %s' % str(sortby))
        ortho.groups.sort(key=key, reverse=(sortby < 0))
    call_function[ENUM_CALLKEY.SORTBY] = call_sortby

    def call_subset_samples(sample_names):
        if sample_names:
            ortho.subset_samples(sample_names, inplace=True)
    call_function[ENUM_CALLKEY.SUBSET_SAMPLES] = call_subset_samples
            
    for call_key, call_input in call_order:
        call_function[call_key](call_input)
    
    return ortho
