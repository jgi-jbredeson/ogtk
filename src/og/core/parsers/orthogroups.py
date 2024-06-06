

from og.constants import (
    _COMMA,
    _COMMENT,
    _EMPTY,
    _EOL,
    _SPACE,
    _TAB,
    range
)
from ..io import open, is_stream

num = len
_CS = _COMMA + _SPACE

def _join_on_comma(l):
    return _EMPTY if l is None else _CS.join(l)


class Orthogroups(object):
    def __init__(self, ):
        self.species = []
        self.ids = []
        self.groups = []


class OrthoFinderOrthogroups(Orthogroups):
    def __init__(self, infile=None, **kwargs):
        Orthogroups.__init__(self)
        if infile is not None:
            if is_stream(infile):
                # io object
                self._parse_orthogroups_file(infile)
            else:
                if 'mode' in kwargs:
                    if 'a' in kwargs['mode'] or \
                       'w' in kwargs['mode']:
                        raise ValueError("%s() constructor is read-only" % (
                            self.__class__.__name__
                        ))
                else:
                    kwargs['mode'] = 'rt'
                self._read_orthogroups_file(open(infile, **kwargs))

                
    def _read_orthogroups_file(self, infile):
        num_fields = -1
        for line in infile:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue
        
            if line.startswith('Orthogroup'):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.species = list(map(str.strip, fields[1:]))
                    
            elif num_fields < 0:
                raise Exception("No header detected in file: %s" % infile)

            else:
                fields = line.split(_TAB)
                group = []
                for i in range(1, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group.append(tuple(map(str.strip, fields[i].split(_COMMA))))
                    else:
                        group.append(tuple())

                self.ids.append(fields[0].strip())
                self.groups.append(group)

        assert num(self.ids) == num(self.groups), \
            "Mismatched number of orthogroups and IDs"

                
    def __iter__(self):
        for i in range(num(self.groups)):
            yield (self.ids[i], self.groups[i])

            
    def format_orthogroups_header(self, species=None, id='Orthogroup'):
        if species is None:
            species = self.species
        return '%s\t%s' % (id, _TAB.join(species))


    def format_orthogroups_record(self, id=None, group=None, index=None):
        if index is not None:
            id = self.ids[index]
            group = self.groups[index]
        return '%s\t%s' % (str(id), _TAB.join(map(_join_on_comma, group)))


    def to_table(self, infile, **kwargs):
        if is_stream(infile):
            stream = infile
            close = False
        else:
            if 'mode' in kwargs:
                if 'r' in kwargs['mode']:
                    raise ValueError("%s.to_table() is write-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'w'
                stream = open(infile, **kwargs)
                close = True

        stream.write(self.format_orthogroups_header() + _EOL)
        for i in range(num(self.groups)):
            stream.write(self.format_orthogroups_record(index=i) + _EOL)
        if close:
            stream.close()
        
