
from og.core.compression import open, is_stream
from og.constants import (
    _COMPRESSION_FLAGS,
    _COMMENT,
    _COMMA,
    _EMPTY,
    _SPACE,
    _TAB,
    range
)
from og.core.orthogroups import (
    Orthogroups,
    Orthogroup,
    Samples,
    Sample, 
)


_CS = _COMMA + _SPACE



def num(items):
    return 0 if None else len(items)



def _join_on_comma(l):
    return _EMPTY if l is None else _CS.join(sorted(l))



class OrthogroupsFormatError(Exception):
    pass



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

                
    def format_header(self, samples=None, id=None, species=None):
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


    def format_record(self, id=None, group=None, index=None):
        if index is None and group is None:
            raise ValueError("group or index keyword argument required")
        elif group is None:
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

    
    def format_record(
            self,
            id=None, cluster=None, prob=None,
            group=None, count=0, index=None
    ):
        if index is None and group is None:
            raise ValueError("group or index keyword argument required")
        elif group is None:
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
        

    def copy(self):
        return self.__class__(
            id=self.id,
            members=self,
            cluster=self.cluster,
            samples=self.samples,
            count=self.count,
            probability=self.probability
        )


        
        
class CountedClusteredOrthogroups(OrthoFinderOrthogroups):
    def __init__(self, file=None, **kwargs):
        kwargs['recordclass'] = kwargs.get('recordclass', CountedClusteredOrthogroup)
        kwargs['prefix'] = kwargs.get('prefix', 'Count\tCluster')
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

                
    def format_header(self, samples=None, id=None, species=None):
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

    
    def format_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is None and group is None:
            raise ValueError("group or index keyword argument required")
        elif group is None:
            group = self.groups[index]
            cluster = group.cluster
            count = group.count
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

        
    def format_header(self, samples=None, id=None, species=None):
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

    
    def format_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is None and group is None:
            raise ValueError("group or index keyword argument required")
        elif group is None:
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



def OrthogroupsFile(filename):
    filehandle = open(filename, 'rt')

    firstline = next(filehandle)
    if firstline.startswith('HOG') or \
       firstline.startswith('Orthogroup'):
        constructor = OrthoFinderOrthogroups
        
    elif firstline.startswith(ClusteredOrthogroups()._prefix):
        constructor = ClusteredOrthogroups

    elif firstline.startswith(CountedClusteredOrthogroups()._prefix):
        if firstline.rstrip().endswith("ManualID\tManualProb\tPostID\tPostProb"):
            constructor = ClusterErrorOrthogroups
        else:
            constructor = CountedClusteredOrthogroups
    else:
        raise OrthogroupsFormatError(
            "Could not determine file format. Header possibly missing"
        )
    
    if filehandle.seekable():
        filehandle.seek(0)
        return constructor(filehandle)
    else:
        ortho = constructor()
        ortho.from_string(firstline + ''.join(filehandle))
        return ortho
