
from og.constants import _EMPTY, _TAB
from og.core.orthogroups.formatters import UpsetFormatter



def upset(
        ortho,
        config=None,
        max_count=0,
        sort_by=0,
        group_by=0,
        force_ascii=False,
        upset_separator=_EMPTY,
        field_separator=_TAB,
        ignore_unplaced=False,
        ignore_unlocalized=None,
        input_reference_names=False,
        map_to_reference_names=False,
):
    if config:
        config.check_samples([s.id for s in ortho.samples])

        if map_to_reference_names:
            ortho.map_loci_to_references(
                config.samples[sample.id],
                is_placed=ignore_unlocalized,
                ignore_unplaced=ignore_unplaced,
                inplace=True
            )

        if map_to_reference_names or \
           input_reference_names:
            ortho.filter_unplaced_references(
                config.samples[sample.id],
                is_placed=ignore_unlocalized,
                ignore_unplaced=ignore_unplaced,
                aggregate_unplaced=True,
                inplace=True
            )
        
    upset = UpsetFormatter(ortho, max_count=max_count)
    upset.group(by=group_by)
    upset.sort(by=sort_by)

    return upset.format(
        ascii=force_ascii,
        field_sep=field_separator,
        upset_sep=upset_separator
    )
