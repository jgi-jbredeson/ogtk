
import sys


from og.constants import (
    _COMMA,
    _COMMENT,
    _EMPTY,
    _EOL,
    _SPACE,
    _TAB,
    range
)
from og.core.compression import open, is_stream

_COMPRESSION_FLAGS = (
    'mode',
    'compresslevel',
    'encoding',
    'errors',
    'newline',
    'compression'
)

_CS = _COMMA + _SPACE



def num(items):
    return 0 if None else len(items)


def _join_on_comma(l):
    return _EMPTY if l is None else _CS.join(l)



class OrthogroupsFormatError(Exception):
    pass



class Sample(object):
    def __init__(self, id, name=None, index=-1):
        self.id = id
        self.name = name
        self.index = index

    def copy(self):
        return self.__class__(self.id, name=self.name, index=self.index)


    
class Samples(list):
    def __init__(self, samples=[]):
        super().__init__()
        self.extend(samples)

        
    def __setitem__(self, index, sample):
        if isinstance(sample, Sample):
            sample = sample.copy()
        else:
            sample = Sample(sample)
        sample.index = index
        super().__setitem__(index, sample)

        
    def append(self, sample):
        if isinstance(sample, Sample):
            sample = sample.copy()
        else:
            sample = Sample(sample)
        sample.index = len(self)
        super().append(sample)


    def extend(self, samples):
        for sample in samples:
            self.append(sample)

            
    def insert(self, index, sample):
        if isinstance(sample, Sample):
            sample = sample.copy()
        else:
            sample = Sample(sample)
        super().insert(index, sample)
        for i, sample in range(index, len(self)):
            self[i].index = i


            
class Members(set):
    def __init__(self, members=(), group=None, sample=None, species=None):
        super().__init__(members)
        self.sample = sample or species
        self.group = group
    


class Orthogroup(list):
    def __init__(self, members=(), id=None, cluster=None, samples=None):
        super().__init__()
        self.cluster = cluster
        self.samples = samples
        self.id = id
        if members:
            self.extend(members)
        elif samples:
            self.extend([()] * len(samples))
        
        
    def __setitem__(self, key, value):
        sample = None
        if self.samples:
            if len(self.samples) <= key:
                raise IndexError(
                    "orthogroup sample count exceeds header sample count"
                )
            sample = self.samples[key] if self.samples else None
        super().__setitem__(
            key,
            Members(value, group=self, sample=sample)
        )

    def __repr__(self):
        return '%s(id=%s, %s)' % (
            self.__class__.__name__,
            self.id,
            super().__repr__()
        )

    def append(self, value):
        sample = None
        if self.samples:
            if len(self.samples) <= len(self):
                raise IndexError(
                    "orthogroup sample count exceeds header sample count"
                )
            sample = self.samples[len(self)]
        super().append(Members(value, group=self, sample=sample))

        
    def extend(self, values):
        for value in values:
            self.append(value)


    @property
    def num_members(self):
        return sum(map(num, self))

    
    @property
    def num_samples(self):
        return sum(map(bool, self))

    
        
class Orthogroups(object):
    def __init__(self, **kwargs):
        self.clear()
        self._prefix  = kwargs.get('prefix', _EMPTY)
        self._factory = kwargs.get('factory', Orthogroup)
        self.filename =	kwargs.get('filename', None)
        self.samples  = kwargs.get('samples', [])
        self.samples  = Samples(self.samples)

        
    def __iter__(self):
        for i in range(num(self.groups)):
            yield self.groups[i]

            
    def format_orthogroups_header(self, samples=None, id=None, species=None):
        raise NotImplementedError('format_orthogroups_header()')

    
    def format_orthogroups_record(self, id=None, group=None, index=None):
        raise NotImplementedError('format_orthogroups_record()')

    
    def to_file(self, file=sys.stdout, **kwargs):
        if is_stream(file):
            stream = file
            close = False
        else:
            if 'r' in kwargs.get('mode', _EMPTY):
                raise ValueError("%s.to_file() is write-only" % (
                    self.__class__.__name__
                ))
            else:
                kwargs['mode'] = 'w'

            kwargs = {
                k:v for k,v in kwargs.items() if k in _COMPRESSION_FLAGS
            }
            stream = open(file, **kwargs)
            close = True

        stream.write(self.format_orthogroups_header() + _EOL)
        for i in range(num(self.groups)):
            stream.write(self.format_orthogroups_record(index=i) + _EOL)

        if close:
            stream.close()


    def from_file(self, file, **kwargs):
        raise NotImplementedError('from_file()')


    def from_string(self, instring):
        raise NotImplementedError('from_string()')


    def clear(self):
        self.samples = []
        self.groups = []
        self.filename = None

        
    def new_group(self, append=False):
        group = self._factory(samples=self.samples)
        if append:
            self.groups.append(group)
        return group

    
    @property
    def species(self):
        return self.samples

    
    @species.setter
    def species(self, species):
        self.samples = samples
        
    to_table = to_file

    

class OrthoFinderOrthogroups(Orthogroups):
    def __init__(self, file=None, **kwargs):
        kwargs['prefix'] = kwargs.get('prefix', 'Orthogroup')
        super().__init__(**kwargs)
        if file is not None:
            self.from_file(file, **kwargs)
        

    def _parse(self, file):
        is_hog = False
        num_fields = -1
        samples_field = 1
        for line in file:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)
            
            if ((num_fields < 0) and 
                ((fields[0] == 'HOG') or 
                 (fields[0] == self._prefix))):
                if fields[0] == 'HOG':
                    samples_field = 3
                    is_hog = True
                num_fields = num(fields)
                self.samples = Samples(map(str.strip, fields[samples_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    str(self.filename)
                ))

            else:
                group = self.new_group(append=True)
                
                fields = line.split(_TAB)
                for i in range(samples_field, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group[i-samples_field].update(
                            map(str.strip, fields[i].split(_COMMA))
                        )
                if is_hog:
                    group.id = fields[0].strip()
                    group.cluster = fields[1].strip()
                else:
                    group.id = fields[0].strip()
                    group.cluster = group.id

                
    def format_orthogroups_header(self, samples=None, id=None, species=None):
        if id is None:
            id = self._prefix
        if samples or species:
            samples = samples or species
        else:
            samples = self.samples
        header = [id]
        for sample in samples:
            if isinstance(sample, Sample):
                sample = sample.id
            header.append(sample)
        return _TAB.join(header)


    def format_orthogroups_record(self, id=None, group=None, index=None):
        if index is not None:
            group = self.groups[index]
            id = group.id
        return '%s\t%s' % (str(id), _TAB.join(map(_join_on_comma, group)))

    
    def from_file(self, file, **kwargs):
        self.clear()
        if is_stream(file):
            # io object
            self.filename = getattr(file, 'name', None)
            self._parse(file)
        else:
            if 'mode' in kwargs:
                if 'a' in kwargs['mode'] or \
                   'w' in kwargs['mode']:
                    raise ValueError("%s() constructor is read-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'rt'

            kwargs = {
                k:v for k,v in kwargs.items() if k in _COMPRESSION_FLAGS
            }
            self.filename = file
            with open(file, **kwargs) as fd:
                self._parse(fd)


    def from_string(self, instring):
        import io
        self.clear()
        self._parse(io.StringIO(instring))

        
    def clear(self):
        super().clear()
        self.filename = None
        
        

class ClusteredOrthogroups(OrthoFinderOrthogroups):
    def __init__(self, file=None, **kwargs):
        kwargs['prefix'] = kwargs.get('prefix','Cluster\tOrthogroup')
        super().__init__(**kwargs)
        if file is not None:
            self.from_file(file, **kwargs)
            
            
    def _parse(self, file):
        cluster_id = -1
        num_fields = -1
        samples_field = 2
        header = self._prefix
        comment = _COMMENT + _SPACE
        for line in file:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(comment):
                continue
            
            if ((num_fields < 0) and
                (line.startswith(header))):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.samples = Samples(map(str.strip, fields[samples_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    getattr(file,'name','<iobuffer>')
                ))

            else:
                group = self.new_group(append=True)
                
                fields = line.split(_TAB)
                for i in range(samples_field, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group[i-samples_field].update(
                            map(str.strip, fields[i].split(_COMMA))
                        )
                group.id = fields[1].strip()                        
                group.cluster = fields[0].strip()

    
    def format_orthogroups_record(
            self,
            id=None, cluster=None, prob=None,
            group=None, count=0, index=None
    ):
        if index is not None:
            group = self.groups[index]
            cluster = group.cluster
            id = group.id
            
        return _TAB.join((
            str(cluster),
            str(id),
            _TAB.join(map(_join_on_comma, group))
        ))
            
        

class CountedClusteredOrthogroup(Orthogroup):
    def __init__(
            self, members=(), id=None, cluster=None,
            samples=None, count=0, probability=-1
    ):
        super().__init__(members, id, cluster, samples)
        self.count = count
        self.probability = probability


        
class CountedClusteredOrthogroups(OrthoFinderOrthogroups):
    def __init__(self, file=None, **kwargs):
        kwargs['prefix'] = kwargs.get('prefix', 'Count\tCluster')
        kwargs['factory'] = kwargs.get('factory', CountedClusteredOrthogroup)
        super().__init__(**kwargs)
        self.grouptag = kwargs.get('grouptag', '##group=')
        if file is not None:
            self.from_file(file, **kwargs)        

            
    def _parse(self, file):
        cluster_id = -1
        num_fields = -1
        samples_field = 2
        group_tag = self.grouptag
        header = self._prefix
        comment = _COMMENT + _SPACE
        for line in file:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(comment):
                continue
            
            if ((num_fields < 0) and
                (line.startswith(header))):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.samples = Samples(map(str.strip, fields[samples_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    getattr(file,'name','<iobuffer>')
                ))

            elif line.startswith(group_tag):
                cluster_id = line[len(group_tag):].strip()
                
            else:
                group = self.new_group(append=True)
                
                fields = line.split(_TAB)
                for i in range(samples_field, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group[i-samples_field].update(
                            map(str.strip, fields[i].split(_COMMA))
                        )
                group.cluster = cluster_id
                group.count = int(fields[0].strip())
                group.id = fields[1].strip()

                
    def format_orthogroups_header(self, samples=None, id=None, species=None):
        if id is None:
            id = self._prefix
        if samples or species:
            samples = samples or species
        else:
            samples = self.samples
        header = [id]
        for sample in samples:
            if isinstance(sample, Sample):
                sample = sample.id
            header.append(sample)
        header.extend(('PostID','PostProb'))
        return _TAB.join(header)

    
    def format_orthogroups_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is not None:
            group = self.groups[index]
            count = group.count
            cluster = group.cluster
            prob = group.probability
            id = group.id
            
        return _TAB.join((
            str(count),
            str(id),
            _TAB.join(map(_join_on_comma, group)),
            str(cluster),
            '%g' % prob
        ))        

        
    
class ClusterErrorOrthogroups(CountedClusteredOrthogroups):
    def __init__(self, file=None, **kwargs):
        self._prefix = kwargs.get('prefix', 'Count\tCluster')        
        self.grouptag = kwargs.get('grouptag', '##group=')
        super().__init__(**kwargs)
        if file is not None:
            self.from_file(file, **kwargs)

            
    def _parse(self, file):
        old_cluster_id = -1
        num_fields = -1
        samples_field = 2
        group_tag = self.grouptag
        header = self._prefix
        comment = _COMMENT + _SPACE
        
        for line in file:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(comment):
                continue
            
            if ((num_fields < 0) and
                (line.startswith(header))):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.samples = Samples(map(str.strip, fields[samples_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    getattr(file,'name','<iobuffer>')
                ))

            elif line.startswith(group_tag):
                old_cluster_id = line[len(group_tag):].strip()
                
            else:
                group = self.new_group(append=True)
                
                fields = line.split(_TAB)
                for i in range(samples_field, num_fields-4):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group[i-samples_field].update(
                            map(str.strip, fields[i].split(_COMMA))
                        )
                old_cluster_id = fields[-4].strip()
                new_cluster_id = fields[-2].strip()
                old_prob = float(fields[-3].strip())
                new_prob = float(fields[-1].strip())

                group.cluster = (old_cluster_id, new_cluster_id)
                group.count = int(fields[0].strip())
                group.id = fields[1].strip()
                group.probability = (old_prob, new_prob)

        
    def format_orthogroups_header(self, samples=None, id=None, species=None):
        if id is None:
            id = self._prefix
        if samples or species:
            samples = samples or species
        else:
            samples = self.samples
        header = [id]
        for sample in samples:
            if isinstance(sample, Sample):
                sample = sample.id
            header.append(sample)
        header.extend(('ManualID','ManualProb','PostID','PostProb'))
        return _TAB.join(header)

    
    def format_orthogroups_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is not None:
            group = self.groups[index]
            cluster = group.cluster
            count = group.count
            prob = group.probability
            id = group.id
        
        return _TAB.join((
            str(count), 
            str(id),
            _TAB.join(map(_join_on_comma, group)),
            str(cluster[0]), '%g' % prob[0],
            str(cluster[1]), '%g' % prob[1],
        ))
