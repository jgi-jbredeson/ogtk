

def update(
        training_ortho, error_ortho,
        min_prob_reassing=-1.0, ignore_groups=()
):
    discard = dict()
    reassign = dict()
    preserve = dict()
    for group in error_ortho.groups:
        if isinstance(group.cluster, (list, tuple)) and len(group.cluster) == 2:
            if group.cluster[1] in ignore_groups:
                continue
        else:
            raise TypeError('orthogroup.cluster must be a 2-tuple')
        
        if isinstance(group.probability, (list, tuple)) and len(group.probability) == 2:
            if group.probability[1] >= min_prob_reassign:
                reassign[group.id] = (group.cluster[1], group.probability[1])
        else:
            raise TypeError('orthogroup.probability must be a 2-tuple')

    for group in training_ortho.groups:
        if group.id in reassign:
            group.cluster = reassign[group.id][0]
            group.probability = reassign[group.id][1]

    return training_ortho
