
from og.constants import (
    _PYTHON_VERSION,
    _COMMENT,
    _EMPTY,
    _TAB
)
from ..io import open, is_stream


if _PYTHON_VERSION < (3,7):
    import collections.OrderedDict as dict

num = len
_STRAND = {'-': -1, '.': 0, '+': +1}
    

class BEDFormatError(Exception):
    pass


class BEDRecord(object):
    def __init__(self, chr, beg, end, name=None, score=0, strand=0):
        self.chr = chr
        self.beg = beg
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand


class BED(dict):
    def __init__(self, infile, **kwargs):
        dict.__init__(self)
        
        if infile is not None:
            self.from_file(infile, **kwargs)


    def _read_file(self, infile):
        line_count = 0
        for line in infile:
            line = line.strip()
            line_count += 1

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)

            if num(fields) < 3:
                raise BEDFormatError(
                    "Expected at least three BED fields, "
                    "line %d" % line_count
                )

            bed = BEDRecord(
                fields[0].strip(),
                int(fields[1]),
                int(fields[2])
            )
            
            if num(fields) > 3:
                bed.name = fields[3].strip()
            if num(fields) > 5:
                bed.score = float(fields[4])
                bed.strand = _STRAND[fields[5].strip()]

            if bed.chr not in self:
                self[bed.chr] = list()

            self[bed.chr].append(bed)


    def from_file(self, infile, **kwargs):
        self.clear()
        if is_stream(infile):
            self._read_file(infile)
        else:
            if 'mode' in kwargs:
                if 'a' in kwargs['mode'] or \
                   'w' in kwargs['mode']:
                    raise ValueError("%s() constructor is read-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'rt'
            with open(infile, **kwargs) as fd:
                self._read_file(fd)


    def from_string(self, instring):
        import io
        self.clear()
        self._read_file(io.StringIO(instring))



class BEDNameMap(BED):
    def _read_file(self, infile):
        line_count = 0
        for line in infile:
            line = line.strip()
            line_count += 1

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)

            if num(fields) < 4:
                raise BEDFormatError(
                    "Expected at least four BED fields, "
                    "line %d" % line_count
                )

            bed = BEDRecord(
                fields[0].strip(),
                int(fields[1]),
                int(fields[2]),
                fields[3].strip()
            )
            
            if num(fields) > 5:
                bed.score = float(fields[4])
                bed.strand = _STRAND[fields[5].strip()]

            if bed.name in self:
                raise KeyError("Duplicate locus name: %s" % bed.name)
            else:
                self[bed.name] = bed

