
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

num = len
_CS = _COMMA + _SPACE


class OrthogroupsFormatError(Exception):
    pass



class Orthogroups(object):
    def __init__(self):
        self.clear()
        self._prefix = ''

        
    def __iter__(self):
        for i in range(num(self.groups)):
            yield (self.ids[i], self.groups[i])

            
    def _join_on_comma(self, l):
        return _EMPTY if l is None else _CS.join(l)
            
            
    def format_orthogroups_header(self, species=None, id=''):
        raise NotImplementedError('format_orthogroups_header()')

    
    def format_orthogroups_record(self, id=None, group=None, index=None):
        raise NotImplementedError('format_orthogroups_record()')

    
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

        stream.write(self.format_orthogroups_header() + _EOL)
        for i in range(num(self.groups)):
            stream.write(self.format_orthogroups_record(index=i) + _EOL)

        if close:
            stream.close()


    def from_file(self, infile, **kwargs):
        raise NotImplementedError('from_file()')


    def from_string(self, instring):
        raise NotImplementedError('from_string()')

    
    def clear(self):
        self.species = []
        self.ids = []
        self.clusters = self.ids
        self.groups = []
        
    to_table = to_file

    

class OrthoFinderOrthogroups(Orthogroups):
    def __init__(self, infile=None, **kwargs):
        super().__init__()
        self._prefix = 'Orthogroup'
        self.is_hog = False
        self.filename = None
        if infile is not None:
            self.from_file(infile, **kwargs)
        

    def _parse(self, infile):
        num_fields = -1
        species_field = 1
        for line in infile:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)
            
            if ((num_fields < 0) and 
                ((fields[0] == 'HOG') or 
                 (fields[0] == self._prefix))):
                if fields[0] == 'HOG':
                    self.clusters = []
                    self.is_hog = True
                    species_field = 3
                num_fields = num(fields)
                self.species = list(map(str.strip, fields[species_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    str(self.filename)
                ))

            else:
                fields = line.split(_TAB)
                group = []
                for i in range(species_field, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group.append(tuple(map(str.strip, fields[i].split(_COMMA))))
                    else:
                        group.append(tuple())
                        
                if self.is_hog:
                    self.clusters.append(fields[1].strip())
                    
                self.ids.append(fields[0].strip())
                self.groups.append(group)

        assert num(self.ids) == num(self.groups), \
            "Mismatched number of orthogroups and IDs"

                
    def format_orthogroups_header(self, species=None, id=None):
        if id is None:
            id = self._prefix
        if species is None:
            species = self.species
        return '%s\t%s' % (id, _TAB.join(species))


    def format_orthogroups_record(self, id=None, group=None, index=None):
        if index is not None:
            id = self.ids[index]
            group = self.groups[index]
        return '%s\t%s' % (str(id), _TAB.join(map(self._join_on_comma, group)))


    def from_file(self, infile, **kwargs):
        self.clear()
        if is_stream(infile):
            # io object
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


    def from_string(self, instring):
        import io
        self.clear()
        self._parse(io.StringIO(instring))

        
    def clear(self):
        super().clear()
        self.filename = None
        self.is_hog = False

        

class ClusteredOrthogroups(OrthoFinderOrthogroups):
    def __init__(self, infile=None, **kwargs):
        super().__init__()
        self._prefix = 'Cluster\tOrthogroup'
        self.clusters = []
        if infile is not None:
            self.from_file(infile, **kwargs)

            
    def _parse(self, infile):
        cluster_id = -1
        num_fields = -1
        species_field = 2
        header = self._prefix
        comment = _COMMENT + _SPACE
        for line in infile:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(comment):
                continue
            
            if ((num_fields < 0) and
                (line.startswith(header))):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.species = list(map(str.strip, fields[species_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    getattr(infile,'name','<iobuffer>')
                ))

            else:
                fields = line.split(_TAB)
                group = []
                for i in range(species_field, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group.append(tuple(map(str.strip, fields[i].split(_COMMA))))
                    else:
                        group.append(tuple())

                self.clusters.append(fields[0].strip())
                self.ids.append(fields[1].strip())
                self.groups.append(group)
                
        assert num(self.clusters) == num(self.ids) == num(self.groups), \
            "Mismatched number of orthogroups, clusters, or IDs"

    
    def format_orthogroups_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is not None:
            id = self.ids[index]
            group = self.groups[index]
            cluster = self.clusters[index]
            
        return _TAB.join((
            str(cluster),
            str(id),
            _TAB.join(map(self._join_on_comma, group))
        ))
            
    def clear(self):
        super().clear()
        self.clusters = []

        
        
class CountedClusteredOrthogroups(OrthoFinderOrthogroups):
    def __init__(self, infile=None, grouptag='##group=', **kwargs):
        super().__init__()
        self._prefix = 'Count\tCluster'
        self.grouptag = grouptag
        self.counts = []
        self.clusters = []
        self.probabilities = []
        if infile is not None:
            self.from_file(infile, **kwargs)
            

    def _parse(self, infile):
        cluster_id = -1
        num_fields = -1
        species_field = 2
        group_tag = self.grouptag
        header = self._prefix
        comment = _COMMENT + _SPACE
        for line in infile:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(comment):
                continue
            
            if ((num_fields < 0) and
                (line.startswith(header))):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.species = list(map(str.strip, fields[species_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    getattr(infile,'name','<iobuffer>')
                ))

            elif line.startswith(group_tag):
                cluster_id = line[len(group_tag):].strip()
                
            else:
                fields = line.split(_TAB)
                group = []
                for i in range(species_field, num_fields):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group.append(tuple(map(str.strip, fields[i].split(_COMMA))))
                    else:
                        group.append(tuple())

                self.clusters.append(cluster_id)
                self.counts.append(int(fields[0].strip()))
                self.ids.append(fields[1].strip())
                self.groups.append(group)
                
        assert num(self.counts) == num(self.clusters) == num(self.ids) == num(self.groups), \
            "Mismatched number of orthogroups, clusters, counts, or IDs"

        self.probabilities = [-1] * num(self.groups)
        
        
    def format_orthogroups_header(self, species=None, id=None):
        if id is None:
            id = self._prefix
        if species is None:
            species = self.species
        return '%s\t%s\tPostID\tPostProb' % (id, _TAB.join(species))

    
    def format_orthogroups_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is not None:
            id = self.ids[index]
            group = self.groups[index]
            count = self.counts[index]
            cluster = self.clusters[index]
            prob = self.probabilities[index]
            
        return _TAB.join((
            str(count),
            str(id),
            _TAB.join(map(self._join_on_comma, group)),
            str(cluster),
            '%g' % prob
        ))        
            
    def clear(self):
        super().clear()
        self.counts = []
        self.clusters = []
        self.probabilities = []


        
class ClusterErrorOrthogroups(CountedClusteredOrthogroups):
    def __init__(self, infile=None, grouptag='##group=', **kwargs):
        super().__init__()
        self._prefix = 'Count\tCluster'
        self.grouptag = grouptag
        self.counts = []
        self.clusters = []
        self.probabilities = []
        self.postclusters = []
        if infile is not None:
            self.from_file(infile, **kwargs)

        
    def _parse(self, infile):
        old_cluster_id = -1
        num_fields = -1
        species_field = 2
        group_tag = self.grouptag
        header = self._prefix
        comment = _COMMENT + _SPACE
        
        for line in infile:
            line = line.lstrip().rstrip('\r\n')

            if line == _EMPTY or \
               line.startswith(comment):
                continue
            
            if ((num_fields < 0) and
                (line.startswith(header))):
                fields = line.split(_TAB)
                num_fields = num(fields)
                self.species = list(map(str.strip, fields[species_field:]))
                
            elif num_fields < 0:
                raise OrthogroupsFormatError("No header detected in: %s" % (
                    getattr(infile,'name','<iobuffer>')
                ))

            elif line.startswith(group_tag):
                old_cluster_id = line[len(group_tag):].strip()
                
            else:
                fields = line.split(_TAB)
                group = []
                for i in range(species_field, num_fields-4):
                    fields[i] = fields[i].strip()
                    if num(fields[i]) > 0:
                        group.append(tuple(map(str.strip, fields[i].split(_COMMA))))
                    else:
                        group.append(tuple())

                old_cluster_id = fields[-4]
                new_cluster_id = fields[-2]
                old_prob = float(fields[-3])
                new_prob = float(fields[-1])
                
                self.clusters.append((old_cluster_id, new_cluster_id))
                self.counts.append(int(fields[0].strip()))
                self.ids.append(fields[1].strip())
                self.groups.append(group)
                self.probabilities.append((old_prob, new_prob))
                
        assert num(self.counts) == num(self.clusters) == num(self.ids) == num(self.groups) == num(self.probabilities), \
            "Mismatched number of orthogroups, clusters, counts, or IDs"


        
    def format_orthogroups_header(self, species=None, id=None):
        if id is None:
            id = self._prefix
        if species is None:
            species = self.species
        return '%s\t%s\tManualID\tManualProb\tPostID\tPostProb' % (id, _TAB.join(species))

    
    def format_orthogroups_record(self, id=None, cluster=None, prob=None, group=None, count=0, index=None):
        if index is not None:
            id = self.ids[index]
            group = self.groups[index]
            count = self.counts[index]
            cluster = self.clusters[index]
            prob = self.probabilities[index]
        
        return _TAB.join((
            str(count), 
            str(id),
            _TAB.join(map(self._join_on_comma, group)),
            str(cluster[0]), '%g' % prob[0],
            str(cluster[1]), '%g' % prob[1],
        ))
