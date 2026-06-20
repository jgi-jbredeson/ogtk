
import sys

from og.core.compression import open, is_stream
from og.core.members import map_loci_to_reference_counts as _map
from og.core.members import filter_unplaced_references as _flt
from og.constants import (
    _COMPRESSION_FLAGS,
    _EMPTY,
    range
)



def num(items):
    return 0 if None else len(items)



class Sample(object):
    def __init__(self, id, name=None, index=-1):
        self.id = id
        self.name = name or id
        self.index = index

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.id == other.id
        return self.id == other

    def __hash__(self):
        return id(self.id)
    
    def __str__(self):
        return self.id

    def __repr__(self):
        return '%s(%s, index=%s)' % (
            self.__class__.__name__,
            repr(self.id),
            repr(self.index)
        )
    
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


    def copy(self):
        return self.__class__(self)


    def get_indices(self):
        return {sample.id: sample.index for sample in self}


            
class Members(set):
    def __init__(self, members=(), group=None, sample=None, species=None):
        super().__init__(members)
        self.sample = sample or species
        self.group = group

    def copy(self):
        return self.__class__(
            member=self,
            group=self.group,
            sample=self.sample,
            species=self.species
        )        


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
            repr(self.id),
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


    def copy(self):
        return self.__class__(
            id=self.id,
            members=self,
            cluster=self.cluster,
            samples=self.samples
        )

    
        
class Orthogroups(object):
    def __init__(self, **kwargs):
        self.clear()
        self._constructor = kwargs.get('recordclass', Orthogroup)
        self._prefix  = kwargs.get('prefix', _EMPTY)
        self.filename =	kwargs.get('filename', None)
        self.samples  = kwargs.get('samples', [])
        self.groups   = kwargs.get('groups', [])

        
    def __iter__(self):
        return iter(self.groups)

            
    def format_header(self, samples=None, id=None, species=None):
        raise NotImplementedError('format_header()')

    
    def format_record(self, id=None, group=None, index=None):
        raise NotImplementedError('format_record()')

    
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

        print(self.format_header(), file=stream)
        for i in range(num(self.groups)):
            print(self.format_record(index=i), file=stream)

        if close:
            stream.close()

            
    to_table = to_file

            
    def from_file(self, file, **kwargs):
        raise NotImplementedError('from_file()')


    def from_string(self, instring):
        raise NotImplementedError('from_string()')


    def clear(self):
        self.samples = []
        self.groups = []
        self.filename = None

        
    def new_group(self, append=False):
        group = self._constructor(samples=self.samples)
        if append:
            self.groups.append(group)
        return group

    
    def new_group_copy(self, other, append=False):
        group = other.copy()
        group.clear()
        group.samples = self.samples
        group.extend([()] * num(self.samples))
        if append:
            self.groups.append(group)
        return group
        
            
    @property
    def samples(self):
        return self._samples

    
    @samples.setter
    def samples(self, samples):
        self._samples = Samples(samples)

        
    @property
    def species(self):
        return self.samples

    
    @species.setter
    def species(self, species):
        self.samples = samples


    def has_sample(self, sample_id):
        return sample_id in self.samples

    
    def has_samples(self, samples):
        return all(map(self.has_sample, samples))

        
    def check_sample(self, sample):
        if not self.has_sample(sample):
            raise KeyError("Sample not found in orthogroups: %s" % str(sample))

        
    def check_samples(self, samples):
        for sample in samples:
            self.check_sample(sample)
            

    def subset_samples(self, samples, min_samples=-1, inplace=False):
        self.check_samples(samples)

        indices = {s.id:s.index for s in self.samples}
        
        copy = self.__class__(samples=samples)
        for group_ in self.groups:
            group = copy.new_group_copy(group_, append=False)
            for sample in copy.samples:
                group[sample.index] = group_[indices[sample.id]]
            if group.num_samples >= min_samples:
                copy.groups.append(group)
        if inplace:
            self.samples = copy.samples
            self.groups = copy.groups
            copy = self
        return copy


    def prepend_sample_prefixes(self, prefix_delim='|', replace_delim=None, inplace=False):
        if replace_delim is None:
            replace_delim = prefix_delim
        copy = self.__class__(samples=self.samples)
        for group_ in self.groups:
            group = copy.new_group_copy(group_, append=True)
            for sample in self.samples:
                prefix = sample.id + prefix_delim
                for member in group_[sample.index]:
                    if not member.startswith(prefix):
                        member = prefix + member.replace(prefix_delim,
                                                         replace_delim)
                    group[sample.index].add(member)
        if inplace:
            self.samples = copy.samples
            self.groups = copy.groups
            copy = self
        return copy


    def remove_sample_prefixes(self, prefix_delim='|', inplace=False):
        copy = self.__class__(samples=self.samples)
        for group_ in self.groups:
            group = copy.new_group_copy(group_, append=True)
            for sample in self.samples:
                prefix = sample.id + prefix_delim
                for member in group_[sample.index]:
                    if member.startswith(prefix):
                        member = member[len(prefix):]
                    group[sample.index].add(member)
        if inplace:
            self.samples = copy.samples
            self.groups = copy.groups
            copy = self
        return copy
    
            
    def map_loci_to_references(self, config, is_placed=None,
                              ignore_unplaced=False, inplace=False):
        """
        Convert an orthogroups table of locus names to reference names.

        The `config` argument is a SampleConfigFile instance.

        Optional keyword arugments `is_placed` takes a function determining
        whether a reference in config is considered placed on an assembled
        molecule. The `ignore_unplaced` keyword indicates whether (and to 
        what degree) to filter unplaced references. A value of `False` or `0`
        indicates no removal, `1` indicates lenient,  `2` indicates strict
        removal.

        The `inplace` argument determines whether the converted orthogroups
        will be stored in the calling Orthogroup object (True) or whether a
        new copy will be returned (False).
        """
        config.check_samples(map(lambda sample: sample.id, self.samples))
        
        copy = self.__class__(samples=self.samples)
        for group_ in self.groups:
            group = group_ \
                if   inplace \
                else copy.new_group_copy(group_, append=True)
            for sample in self.samples:
                group[sample.index] = _map(
                    group_[sample.index],
                    config.samples[sample.id],
                    is_placed,
                    ignore_unplaced
                )
        return self if inplace else copy


    def filter_unplaced_references(self, config,
                                  is_placed=False, ignore_unplaced=False,
                                  aggregate_unplaced=True, inplace=False):
        """
        Filter an orthogroups table of reference names.

        The `config` argument is a SampleConfigFile instance.

        Optional keyword arugments `is_placed` takes a function determining
        whether a reference in config is considered placed on an assembled
        molecule. The `ignore_unplaced` keyword indicates whether (and to 
        what degree) to filter unplaced references. A value of `False` or `0`
        indicates no removal, `1` indicates lenient,  `2` indicates strict
        removal. `aggregate_unplaced` specifies whether an unplaced reference
        name will be mapped to the `unplaced_id` value specified for that 
        sample in config.

        The `inplace` argument determines whether the filtered orthogroups
        will be stored in the calling Orthogroup object (True) or whether a
        new copy will be returned (False).
        """        
        config.check_samples(map(lambda sample: sample.id, self.samples))

        copy = self.__class__(samples=self.samples)
        for group_ in self.groups:
            group = group_ \
                if   inplace \
                else copy.new_group_copy(group_, append=True)            
            for sample in self.samples:
                group[sample.index] = _flt(
                    group_[sample.index],
                    config.samples[sample.id],
                    is_placed,
                    ignore_unplaced,
                    aggregate_unplaced
                )
        return self if inplace else copy
    

    
        
        

