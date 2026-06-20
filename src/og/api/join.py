
import sys

from og.core.generic import Enumerative
from og.core.members import map_members_to_group_index
from og.core.orthogroups.parsers import OrthoFinderOrthogroups


ENUM_JOIN = Enumerative(INNER=0, LEFT=1, RIGHT=2, FULL=3)
DISJOINT = 'DISJOINT'
NULL = -1


num = len


def format_id(left_id, rght_id, sep=':', reverse=False):
    if reverse:
        return str(rght_id) + sep + str(left_id)
    else:
        return str(left_id) + sep + str(rght_id)

    

def _get_inner_join_indices(
        left_ortho,
        rght_ortho,
        method=ENUM_JOIN.INNER,
        max_left=sys.maxsize,
        max_right=sys.maxsize,
        min_samples=1,
        max_samples=sys.maxsize,
        min_members=1,
        max_members=sys.maxsize
):
    left_samples, innr_samples, rght_samples = \
        _sample_intersection(left_ortho, rght_ortho, method=method)
    
    left_sample_indices = left_ortho.samples.get_indices()
    rght_sample_indices = rght_ortho.samples.get_indices()

    if method == ENUM_JOIN.RIGHT:
        innr_sample_indices = rght_sample_indices
        innr_ortho = rght_ortho
    else:
        innr_sample_indices = left_sample_indices
        innr_ortho = left_ortho

    innr_member_indices = \
        map_members_to_group_index(innr_ortho, samples=innr_samples)
    left_member_indices = \
        map_members_to_group_index(left_ortho)
    rght_member_indices = \
        map_members_to_group_index(rght_ortho)

    
    innr_counts = dict()
    innr_indices = list()
    left_mapped_counts = [0] * num(left_ortho.groups)
    rght_mapped_counts = [0] * num(rght_ortho.groups)
    for member in innr_member_indices:
        if member in left_member_indices and \
           member in rght_member_indices:
            l = left_member_indices[member]
            r = rght_member_indices[member]
            
            if (l,r) in innr_counts:
                innr_counts[(l,r)][1] += 1
            else:
                innr_counts[(l,r)] = [0, 1]

    for l,r in innr_counts:
        innr_counts[(l,r)][0] = sum(bool(
            left_ortho.groups[l][left_sample_indices[s.id]] & \
            rght_ortho.groups[r][rght_sample_indices[s.id]]
        ) for s in innr_samples)

    for l,r in sorted(innr_counts, key=innr_counts.get, reverse=True):
        if ((left_mapped_counts[l] < max_left) and \
            (rght_mapped_counts[r] < max_right)):

            num_samples = innr_counts[(l,r)][0]
            num_members = innr_counts[(l,r)][1]
        
            if ((min_samples <= num_samples <= max_samples) and \
                (min_members <= num_members <= max_members)):
                left_mapped_counts[l] += 1
                rght_mapped_counts[r] += 1
                innr_indices.append((l,r))

    return innr_indices



def _fill_full_outer_indices(inner, Nl, Nr, null=NULL):
    nl = sorted(map(lambda i: i[0], inner))
    nr = sorted(map(lambda i: i[1], inner))
    Ni = num(inner)
    ii = 0
    il = 0
    ir = 0
    jl = 0
    jr = 0
    indices = list()
    inner = sorted(inner)
    while ii < Ni and il < Nl and ir < Nr:
        if il < nl[jl] and ir < nr[jr]:
            if il < ir:
                indices.append((il, null))
                il += 1
            elif il > ir:
                indices.append((null, ir))
                ir += 1
            else:
                indices.append((il, null))
                indices.append((null, ir))
                il += 1
                ir += 1
        elif il < nl[jl]:
            indices.append((il, null))
            il += 1
        elif ir < nr[jr]:
            indices.append((null, ir))
            ir += 1
        else:
            indices.append(inner[ii])
            ii += 1
            if il == nl[jl]:
                il += 1
                jl += 1
            if ir == nr[jr]:
                ir += 1
                jr += 1
            
    while il < Nl or ir < Nr:
        if il < Nl and ir < Nr:
            if il < ir:
                indices.append((il, null))
                il += 1
            elif il > ir:
                indices.append((null, ir))
                ir += 1
        elif il < Nl:
            indices.append((il, null))
            il += 1
        elif ir < Nr:
            indices.append((null, ir))
            ir += 1

    return indices



def _fill_outer_indices(inner, N, side=1, null=NULL):
    """Assumes `inner` is sorted ascending by the values of `side`"""
    Ni = len(inner)
    x = side - 1
    ii = 0
    ix = 0
    indices = list()
    inner = sorted(inner, key=lambda i: i[x])
    while ii < Ni and ix < N:
        if ix < inner[ii][x]:
            indices.append((null, ix) if x else (ix, null))
            ix += 1
        else:
            indices.append(inner[ii])
            if ix == inner[ii][x]:
                ix += 1
            ii += 1
    while ix < N:
        indices.append((null, ix) if x else (ix, null))
        ix += 1
        
    return indices        
        


def _backfill_ortho(join_ortho, innr_samples, join_member_indices):
    innr_member_indices = \
        map_members_to_group_index(join_ortho, samples=innr_samples)
    for member in join_member_indices:
        if member not in innr_member_indices:        
            if num(join_member_indices[member]) == 1:
                g,i = join_member_indices[member].pop()
                join_ortho.groups[g][i].add(member)
            else:
                # new group
                pass
            


def _fill_ortho(
        join_indices,
        left_ortho, rght_ortho,
        method=ENUM_JOIN.INNER,
        strict_intersection=False,
        reverse_id_order=False
):
    left_samples, innr_samples, rght_samples = \
        _sample_intersection(left_ortho, rght_ortho)

    join_ortho = \
        OrthoFinderOrthogroups(
            samples=(
                left_samples + \
                innr_samples + \
                rght_samples
            )
        )
    left_sample_indices = left_ortho.samples.get_indices()
    rght_sample_indices = rght_ortho.samples.get_indices()
    join_sample_indices = join_ortho.samples.get_indices()

    join_index = -1
    left_backfill_members = dict()
    rght_backfill_members = dict()
    for l,r in join_indices:
        join_group = join_ortho.new_group(append=True)
        join_index += 1
        
        left_ogid = DISJOINT
        rght_ogid = DISJOINT
        if l != NULL and r != NULL:
            for sample in innr_samples:
                il = left_sample_indices[sample.id]
                ir = rght_sample_indices[sample.id]
                ij = join_sample_indices[sample.id]

                join_group[ij] = \
                    left_ortho.groups[l][il] & rght_ortho.groups[r][ir]
                
                if not strict_intersection:
                    for member in (left_ortho.groups[l][il] - join_group[ij]):
                        if member in left_backfill_members:
                            left_backfill_members[member].append(
                                (join_index, ij)
                            )
                        else:
                            left_backfill_members[member] = [
                                (join_index, ij)
                            ]
                    for member in (rght_ortho.groups[r][ir] - join_group[ij]):
                        if member in rght_backfill_members:
                            rght_backfill_members[member].append(
                                (join_index, ij)
                            )
                        else:
                            rght_backfill_members[member] = [
                                (join_index, ij)
                            ]
        elif l != NULL:
            for sample in innr_samples:
                il = left_sample_indices[sample.id]
                ij = join_sample_indices[sample.id]
                join_group[ij] = left_ortho.groups[l][il]

        elif r != NULL:
            for sample in innr_samples:
                ir = rght_sample_indices[sample.id]
                ij = join_sample_indices[sample.id]
                join_group[ij] = rght_ortho.groups[r][ir]

        if l != NULL:
            left_ogid = left_ortho.groups[l].id
            for sample in left_samples:
                il = left_sample_indices[sample.id]
                ij = join_sample_indices[sample.id]
                join_group[ij] = left_ortho.groups[l][il]

        if r != NULL:
            rght_ogid = rght_ortho.groups[r].id
            for sample in rght_samples:
                ir = rght_sample_indices[sample.id]
                ij = join_sample_indices[sample.id]
                join_group[ij] = rght_ortho.groups[r][ir]

        join_group.id = \
            format_id(left_ogid, rght_ogid, reverse=reverse_id_order)

        
    if not strict_intersection:
        if method & ENUM_JOIN.LEFT:
            _backfill_ortho(join_ortho, innr_samples, left_backfill_members)
        if method & ENUM_JOIN.RIGHT:
            _backfill_ortho(join_ortho, innr_samples, rght_backfill_members)
        
    return join_ortho
        


def _sample_intersection(left_ortho, rght_ortho, method=ENUM_JOIN.INNER):
    innr_ortho = rght_ortho if method == ENUM_JOIN.RIGHT else left_ortho

    left_sample_indices = left_ortho.samples.get_indices()
    rght_sample_indices = rght_ortho.samples.get_indices()

    innr_samples = set(left_sample_indices) & set(rght_sample_indices)
    left_samples = set(left_sample_indices) - innr_samples
    left_samples = list(
        filter((lambda s: s.id in left_samples), left_ortho.samples)
    )
    rght_samples = set(rght_sample_indices) - innr_samples
    rght_samples = list(
        filter((lambda s: s.id in rght_samples), rght_ortho.samples)
    )
    innr_samples = list(
        filter((lambda s: s.id in innr_samples), innr_ortho.samples)
    )
    return (left_samples, innr_samples, rght_samples)



def join(
        left_ortho,
        rght_ortho,
        method=ENUM_JOIN.INNER,
        max_intersecting_left=sys.maxsize,
        max_intersecting_right=sys.maxsize,
        min_intersecting_members=1,
        max_intersecting_members=sys.maxsize,
        min_intersecting_samples=1,
        max_intersecting_samples=sys.maxsize,
        output_intersecting_strict=False,
        reverse_id_order=False
):
    join_indices = \
        _get_inner_join_indices(
            left_ortho,
            rght_ortho,
            max_left=max_intersecting_left,
            max_right=max_intersecting_right,
            min_samples=min_intersecting_samples,
            max_samples=max_intersecting_samples,
            min_members=min_intersecting_members,
            max_members=max_intersecting_members
        )
    
    if method == ENUM_JOIN.FULL:
        join_indices = \
            _fill_full_outer_indices(
                inner=join_indices,
                Nl=num(left_ortho.groups),
                Nr=num(rght_ortho.groups),
                null=NULL
            )
    elif method == ENUM_JOIN.LEFT:
        join_indices = \
            _fill_outer_indices(
                inner=join_indices,
                N=num(left_ortho.groups),
                side=ENUM_JOIN.LEFT,
                null=NULL
            )
    elif method == ENUM_JOIN.RIGHT:
        join_indices = \
            _fill_outer_indices(
                inner=join_indices,
                N=num(rght_ortho.groups),
                side=ENUM_JOIN.RIGHT,
                null=NULL
            )
    else:
        join_indices.sort()

    return _fill_ortho(
        join_indices,
        left_ortho, rght_ortho, method=method,
        strict_intersection=output_intersecting_strict,
        reverse_id_order=reverse_id_order
    )
