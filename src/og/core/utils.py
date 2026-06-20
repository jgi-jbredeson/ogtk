
from og.constants import dict


def sign(value):
    return value / abs(value)


def count_items(iterable):
    count = dict()
    for item in iterable:
        try:
            count[item] += 1
        except KeyError:
            count[item] = 1
    return count


def index_list(iterable, start=0):
    return dict(zip(iterable, range(start,start+len(iterable))))


def uniq_list(iterable):
    seen = set()
    uniq = list()
    for item in iterable:
        if item not in seen:
            uniq.append(item)
            seen.add(item)
    return uniq


def isiterable(item):
    return hasattr(item, '__iter__') or hasattr(item, '__next__')
