
import sys

from og.core.orthogroups.parsers import OrthoFinderOrthogroups


num = len


def passes_tree_filter(tree, counts, sample_indices):
    passes = True
    # anc = ancestor, dsc = descendant
    for anc_index, dsc_index in tree.get_edges(indices=True, reverse=True):
        dsc = tree.nodes[dsc_index]
        if dsc.length.minimum <= counts[dsc_index] <= dsc.length.maximum:
            if anc_index is not None:
                counts[anc_index] += int(counts[dsc_index] > 0) \
                    if   (dsc.id in sample_indices) \
                    else counts[dsc_index]
        else:
            passes = False
    return passes



def count_tree_members(tree, group, sample_indices):
    counts = [0] * num(tree.nodes)
    for node, index in zip(tree.terminal_nodes, tree.terminal_indices):
        counts[index] = num(group[sample_indices[node.id]])
    return counts



def count_group_members(group):
    return [num(m) for m in group]



def filter(
        ortho, config=None,
        invert=False,
        min_samples=1,
        max_samples=sys.maxsize,
        min_members=1,
        max_members=sys.maxsize,
        ploidy_tree_filter=False,
        ignore_unplaced=False,
        ignore_unlocalized=False,
        input_reference_names=False,
        map_to_reference_names=False,
        output_reference_names=False
):
    tree = None
    if ploidy_tree_filter:
        tree = config.tree['ploidy']
        ortho.check_samples([n.id for n in tree.terminal_nodes])

    if map_to_reference_names:
        mapped = ortho.map_loci_to_references(
            config,
            is_placed=ignore_unlocalized,
            ignore_unplaced=ignore_unplaced
        )
        if output_reference_names:
            ortho = mapped
    else:
        mapped = ortho

    if map_to_reference_names or input_reference_names:
        mapped.filter_unplaced_references(
            config,
            is_placed=ignore_unlocalized,
            ignore_unplaced=ignore_unplaced,
            aggregate_unplaced=True,
            inplace=True
        )

    flt_ortho = OrthoFinderOrthogroups(samples=mapped.samples)
    sample_indices = flt_ortho.samples.get_indices()
    for i, group in enumerate(mapped.groups):
        if ploidy_tree_filter:
            counts = count_tree_members(tree, group, sample_indices)
        else:
            counts = count_group_members(group)

        passes = True
        if ploidy_tree_filter:
            passes = passes_tree_filter(tree, counts, sample_indices)
            
        if not (min_members <= sum(counts) <= max_members):
            passes = False
        if not (min_samples <= sum(map(bool, counts)) <= max_samples):
            passes = False

        if invert:
            passes = not passes
                
        if passes:
            flt_ortho.groups.append(ortho.groups[i])

    return flt_ortho
