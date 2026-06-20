
num = len


from og.core.orthogroups.formatters import OrthogroupCountsTable
    

def count(
        ortho,
        config=None,
        count_samples=False,
        input_reference_names=False,
        map_to_reference_names=False,
):
    if config:
        if map_to_reference_names:
            ortho.map_loci_to_references(
                config,
                is_placed=ignore_unlocalized,
                ignore_unplaced=ignore_unplaced,
                inplace=True
            )

        if map_to_reference_names or input_reference_names:
            ortho.filter_unplaced_references(
                config,
                is_placed=ignore_unlocalized,
                ignore_unplaced=ignore_unplaced,
                aggregate_unplaced=True,
                inplace=True
            )

    return OrthogroupCountsTable(ortho, count_samples=count_samples)
        
