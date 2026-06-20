
num = len

from og.core.members import map_members_to_group_index


def add_singletons(
        ortho, config,
        ignore_unlocalized=None,
        output_reference_names=False,
        ogid_formatter='SGL{0:06d}'.format
):
    member_indices = map_members_to_group_index(ortho)

    singleton_count = 1
    for sample in ortho.samples:
        for locus_name in config.samples[sample.id].loci:
            if locus_name not in member_indices:
                group = ortho.new_group(append=True)
                group.id = ogid_formatter(singleton_count)
                group[sample.index] = (locus_name,)
                singleton_count += 1

    if output_reference_names:
        ortho.map_loci_to_references(
            config.samples[sample.id],
            ignore_unlocalized=ignore_unlocalized,
            ignore_unplaced=False
        )

    return ortho
