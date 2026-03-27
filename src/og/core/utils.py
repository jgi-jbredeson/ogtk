
from og.constants import dict

def count_items(items):
    count = dict()
    for item in items:
        if item in count:
            count[item] += 1
        else:
            count[item] = 1
    return count


def index_list(list_, start=0):
    return dict(zip(list_, range(start,start+len(list_))))


def sorted_orthogroups_by_ids(ortho, key=None, reverse=False):
    if key is None:
        key = lambda item:item
    index = index_list(ortho.ids)
    sorted_ortho = ortho.__class__()
    sorted_ortho.species = ortho.species
    for gid in sorted(ortho.ids, key=key, reverse=reverse):
        sorted_ortho.ids.append(gid)
        sorted_ortho.groups.append(ortho.groups[index[gid]])
    return sorted_ortho


