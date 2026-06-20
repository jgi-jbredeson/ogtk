

from math import inf, isinf
from string import (
    ascii_uppercase,
    ascii_lowercase,
    digits,
    punctuation
)
from collections import deque
from og.core.utils import sign
from og.core.counters import CountsTable
from og.constants import (
    _EMPTY,
    _EOL,
    _TAB
)


GROUP_MEMBERSHIP = 0x1
GROUP_MULTIPLES = 0x2
SORT_CLUSTERS = 0x1
SORT_MEMBERS = 0x2


_LABELS = \
    ascii_uppercase + \
    ascii_lowercase + \
    digits + punctuation

_MAXCHAR = '\u25CF'
_MINCHAR = '\u2015'

_UNICODE_NUMMAP = {
    1 :'\u2460',
    2 :'\u2461',
    3 :'\u2462',
    4 :'\u2463',
    5 :'\u2464',
    6 :'\u2465',
    7 :'\u2466',
    8 :'\u2467',
    9 :'\u2468',
    # 10:'\u2469',
    # 11:'\u246A',
    # 12:'\u246B',
    # 13:'\u246C',
    # 14:'\u246D',
    # 15:'\u246E',
    # 16:'\u246F',
    # 17:'\u2470',
    # 18:'\u2471',
    # 19:'\u2472',
    # 20:'\u2473',
}


num = len


def get_sample_obj(ortho, sample_id):
    if isinstance(sample_id, (bytes, str)):
        for sample in ortho.samples:
            if sample.id == sample_id:
                return sample
    return sample_id



class OrthogroupCountsTable(CountsTable):
    def __init__(self, ortho, count_samples=False):
        super().__init__()
        self._ortho = ortho
        self._init_from_ortho()

    def _init_from_ortho(self, count_samples=False):
        self.samples = self._ortho.samples
        for group in self._ortho.groups:
            record = self.new_record(append=True)
            record.id = group.id
            for sample in self._ortho.samples:
                record[sample.index] = int(bool(group[sample.index])) \
                    if   count_samples \
                    else len(group[sample.index])



class OrthoVennFormatter(object):
    def __init__(self, ortho):
        self._ortho = \
            ortho.__class__(
                samples=ortho.samples,
                groups=list(ortho.groups)
            )
        self._ortho.groups.sort(
            key=lambda g: (g.num_samples, g.num_members),
            reverse=True
        )
        
    def format_header(self, samples=None, id=None, species=None):
        return ''
        
    def format_record(self, group=None, index=None):
        if index is None and group is None:
            raise ValueError("group or index keyword argument required")
        elif group is None:
            group = self._ortho.groups[index]
            
        return _TAB.join([
            str(member) for i in range(len(group)) for member in group[i]
        ])
    
    def __iter__(self):
        for group in self._ortho.groups:
            yield self.format_record(group=group)



class MCScanAnchorsFormatter(object):
    def __init__(self, ortho, sample1, sample2):
        self._ortho = ortho
        self._sample1 = get_sample_obj(ortho, sample1)
        self._sample2 = get_sample_obj(ortho, sample2)

    def format_header(self, samples=None, id=None, species=None):
        return ''

    def format_record(self, group=None, index=None):
        if index is None and group is None:
            raise ValueError("group or index keyword argument required")
        elif group is None:
            group = self._ortho.groups[index]

        if self._sample1 in ("HOG","Orthogroup"):
            members_i = (group.id,)
        else:
            members_i = group[self._sample1.index]
        if self._sample2 in ("HOG","Orthogroup"):
            members_j = (group.id,)
        else:
            members_j = group[self._sample2.index]

        for member_i in members_i:
            for member_j in members_j:
                yield _TAB.join((str(member_i), str(member_j)))

    def __iter__(self):
        for group in self._ortho.groups:
            for record in self.format_record(group=group):
                yield record



def _fmtchar(x):
    try:
        return _UNICODE_NUMMAP[x]
    except KeyError:
        return str(x)



def _fmtcharunicode(x, _min, _max):
    return _min if x < 1 else _max if isinf(x) else _fmtchar(x)



def _fmtcharascii(x, _min, _max):
    return _min if x < 1 else _max if isinf(x) else str(x)



def _flattened(item):
    flattened = []
    queue = deque()
    queue.append(item)
    while queue:
        item = queue.popleft()
        if isinstance(item, list):
            queue.extend(item)
        else:
            flattened.append(item)
    return flattened



def _recursive_sort(item, key=None, reverse=False):
    queue = deque()
    queue.append(item)
    while queue:
        item = queue.popleft()
        if isinstance(item, UpsetGroup):
            item.sort(key=key, reverse=reverse)
            for subitem in item:
                if isinstance(subitem, UpsetGroup):
                    queue.append(subitem)

                    

class UpsetSample(object):
    def __init__(self, id, label=None):
        self.id = id
        self.label = label

        

class UpsetRecord(object):
    def __init__(self, pattern, num_groups=0, num_members=0):
        self.pattern = pattern
        self.num_groups = num_groups
        self.num_members = num_members

    def __hash__(self):
        return id(self.pattern)
        
    def __repr__(self):
        return "%s(%s, %s, %s)" % (
            self.__class__.__name__,
            str(self.pattern),
            str(self.num_groups),
            str(self.num_members)
        )

    def __gt__(self, other):
        return self.pattern > other.pattern

    def __eq__(self, other):
        return self.pattern == other.pattern
    
    

class UpsetGroup(list):
    def __init__(self, iterable=(), num_groups=0, num_members=0):
        super().__init__(iterable)
        self.num_groups = num_groups
        self.num_members = num_members

        

class UpsetFormatter(object):
    def __init__(self, ortho, max_count=0, ascii=False, min_char=_MINCHAR, max_char=_MAXCHAR):
        self.ascii = ascii
        self.max_count = max_count
        self.min_char = min_char
        self.max_char = max_char
        self._ortho = ortho
        self._samples = []
        self._records = []
        self._init_from_ortho()

        
    def _init_from_ortho(self, count_samples=False):
        num_samples = num(self._ortho.samples)
        if num_samples > num(_LABELS):
            raise IndexError("Too many samples")
        
        for i, sample in enumerate(self._ortho.samples):
            self._samples.append(UpsetSample(sample.id, label=_LABELS[i]))

        records = dict()
        for record in OrthogroupCountsTable(self._ortho):
            total = record.total
            for sample in self._ortho.samples:
                record[sample.index] = inf \
                    if  record[sample.index] > self.max_count \
                    else record[sample.index]
            pattern = tuple(record)

            if pattern not in records:
                records[pattern] = UpsetRecord(pattern)
            records[pattern].num_members += total
            records[pattern].num_groups += 1

        self._records = UpsetGroup([UpsetGroup(records.values())])
                

    def sort(self, by=0):
        if abs(by) & SORT_CLUSTERS:
            key = lambda r: (r.num_groups, r.num_members)
        elif abs(by) & SORT_MEMBERS:
            key = lambda r: (r.num_members, r.num_groups)
        elif by:
            raise AssertionError("Invalid sort value: %s" % str(by))
        else:
            key = None

        _recursive_sort(self._records, key=key, reverse=(by < 0))
                

    def group(self, by=0):
        if by & GROUP_MEMBERSHIP:
            if by & GROUP_MULTIPLES:
                key = lambda r: sum(map(isinf, r.pattern))
            else:
                key = lambda r: tuple(sorted(r.pattern))
        else:
            key = lambda r: 0

        prev = None
        group = None
        groups = UpsetGroup()
        for record in sorted(_flattened(self._records), key=key):
            curr = key(record)
            
            if prev != curr:
                groups.append(UpsetGroup())
                group = groups[-1]
            group.append(record)
            group.num_members += record.num_members
            group.num_groups += record.num_groups
            prev = curr
        
        self._records = groups
            
            
    def format(self, ascii=False, min_char=None, max_char=None, upset_sep=_EMPTY, field_sep=_TAB):
        min_char = min_char or self.min_char
        max_char = max_char or self.max_char
        ascii = ascii or self.ascii

        if ascii or self.max_count > num(_UNICODE_NUMMAP):
            fmtchar = lambda x: _fmtcharascii(x, min_char, max_char)
        else:
            fmtchar = lambda x: _fmtcharunicode(x, min_char, max_char)
        
        lines = []
        for sample in self._samples:
            lines.append(field_sep.join([
                "##SAMPLE","{0.label}:{0.id}".format(sample)
            ]))
            
        lines.append("##UPSET")
        lines.append(field_sep.join([
            upset_sep.join((sample.label for sample in self._samples)),
            'Clusters',
            'Proteins'
        ]))
        
        for group in self._records:
            for record in group:
                lines.append(field_sep.join([
                    upset_sep.join(map(fmtchar, record.pattern)),
                    str(record.num_groups),
                    str(record.num_members)
                ]))
            lines.append(field_sep.join([
                "##total",
                str(group.num_groups),
                str(group.num_members)
            ]))

        return _EOL.join(lines)


    def __str__(self):
        return self.format()


    
