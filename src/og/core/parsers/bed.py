
import sys
from og.constants import (
    _PYTHON_VERSION,
    _COMMENT,
    _EMPTY,
    _DOT,
    _EOL,
    _TAB
)
from og.core.compression import open, is_stream


if _PYTHON_VERSION < (3,7):
    from collections import OrderedDict as dict

num = len
_STRAND_TO_INT = {'-': -1, '.': 0, '+': +1}
_STRAND_TO_STR = ('.','+','-')    


class BEDFormatError(Exception):
    pass



class BEDRecord(object):
    def __init__(self, chr, beg, end):
        self.chr = chr
        self.beg = beg
        self.end = end

    def __str__(self):
        return _TAB.join(map(str, (
            self.chr,
            self.beg,
            self.end
        )))

    

class BED6Record(BEDRecord):
    def __init__(self, chr, beg, end, name=None, score=0, strand=0):
        BEDRecord.__init__(self, chr, beg, end)
        self.name = name
        self.score = score
        self.strand = strand

    def __str__(self):
        return _TAB.join(map(str, (
            BEDRecord.__str__(self),
            _DOT if self.name is None else self.name,
            self.score,
            _STRAND_TO_STR[self.strand]
        )))


    
class BEDFile(dict):
    def __init__(self, infile, recordclass=BED6Record, **kwargs):
        dict.__init__(self)
        self.filename = None
        self.max_fields = None
        self.recordclass = recordclass
        if infile is not None:
            self.from_file(infile, **kwargs)


    def _parse(self, infile):
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
            if self.max_fields is None:
                self.max_fields = num(fields)
            if num(fields) != self.max_fields:
                raise BEDFormatError(
                    "Expected %d BED fields, line %d" % (
                        self.max_fields, line_count
                    )
                )
            
            bed = self.recordclass(
                fields[0].strip(),
                int(fields[1]),
                int(fields[2])
            )
            
            if num(fields) > 3:
                bed.name = fields[3].strip()
            if num(fields) > 5:
                bed.score = float(fields[4])
                bed.strand = _STRAND_TO_INT[fields[5].strip()]

            if bed.chr not in self:
                self[bed.chr] = list()

            self[bed.chr].append(bed)


    def from_file(self, infile, **kwargs):
        self.clear()
        if is_stream(infile):
            self.filename = getattr(infile, 'name', None)
            self._parse(infile)
        else:
            if 'mode' in kwargs:
                if 'a' in kwargs['mode'] or \
                   'w' in kwargs['mode']:
                    raise ValueError("%s() constructor is read-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'rt'

            self.filename = infile
            with open(infile, **kwargs) as fd:
                self._parse(fd)


    def to_file(self, file=sys.stdout, **kwargs):
        if is_stream(file):
            stream = file
            close = False
        else:
            if 'mode' in kwargs:
                if 'r' in kwargs['mode']:
                    raise ValueError("%s.to_file() is write-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'w'
                stream = open(file, **kwargs)
                close = True

        for chr in self:
            for record in self[chr]:
                stream.write(str(record) + _EOL)
                
        if close:
            stream.close()
                

                
    def from_string(self, instring):
        import io
        self.clear()
        self._parse(io.StringIO(instring))


    def clear(self):
        dict.clear(self)
        self.filename = None


        
class BEDNameMapFile(BEDFile):
    def _parse(self, infile):
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

            bed = self.recordclass(
                fields[0].strip(),
                int(fields[1]),
                int(fields[2]),
                fields[3].strip()
            )
            
            if num(fields) > 5:
                bed.score = float(fields[4])
                bed.strand = _STRAND_TO_INT[fields[5].strip()]

            if bed.name in self:
                raise KeyError("Duplicate locus name: %s" % bed.name)
            else:
                self[bed.name] = bed

