
import yaml  # requires PyYAML specifically: https://pypi.org/project/PyYAML

from math import inf as _POS_INF
from og.core.bed import BEDNameMapFile
from og.core.newick import IntervalNewickTree
from og.core.assembly_report import AssemblyReportFile
from og.core.compressio import open, is_stream
from og.constants import _PYTHON_VERSION, dict



_NEG_INF = -1.0 * _POS_INF
    

def to_bool(value):
    tmpvalue = value.lower()
    if tmpvalue == 'on' or \
       tmpvalue == 'y' or tmpvalue == 'yes' or \
       tmpvalue == 't' or tmpvalue == 'true':
        return True
    elif tmpvalue == 'off' or \
         tmpvalue == 'n' or tmpvalue == 'no' or \
         tmpvalue == 'f' or tmpvalue == 'false':
        return False
    try:
        tmpvalue = float(tmpvalue)
    except ValueError:
        raise TypeError('Invalid boolean value: %s' % value) from None
    return bool(tmpvalue)


class _SampleRecord(object):
    def __init__(self, sample=None, loci=None, references=None, unplaced_id=None, is_outgroup=False, species=None):
        self.sample = sample or species
        self.loci = loci
        self.references = references
        self.unplaced_id = unplaced_id
        self.is_outgroup = is_outgroup
        

class SampleConfigFile(object):
    def __init__(self, infile, load_files=False, map_assigned_molecule=False):
        self.clear()
        self.filename = None
        self.map_assigned_molecule=map_assigned_molecule
        if infile is not None:
            self.from_file(infile)
        if load_files:
            self.load_files()

            
    def _parse(self, infile):
        yamldata = yaml.safe_load(infile)

        if 'tree' in yamldata:
            self.tree = dict()
            for key in yamldata['tree']:
                self.tree[key] = IntervalNewickTree(yamldata['tree'][key], default_length=(_NEG_INF, _POS_INF))
        else:
            raise ValueError('`tree` key not defined')
        
        for section in yamldata:
            if section.lower() == 'tree':
                continue
            else:  # a sample ID
                self.samples[section] = _SampleRecord(sample=section)
                if 'loci' in yamldata[section]:
                    self.samples[section].loci = yamldata[section]['loci']
                else:
                    raise ValueError('`loci` key not defined for %s' % section)
                
                if 'references' in yamldata[section]:
                    self.samples[section].references = yamldata[section]['references']
                else:
                    raise ValueError('`references` key not defined for %s' % section)
                
                if 'unplaced_id' in yamldata[section]:
                    self.samples[section].unplaced_id = yamldata[section]['unplaced_id']
                else:
                    raise ValueError('`unplaced_id` key not defined for %s' % section)
                
                if 'outgroup' in yamldata[section]:
                    self.samples[section].is_outgroup = to_bool(str(yamldata[section]['outgroup']))

                if self.samples[section].is_outgroup:
                    self._outgroup.append(section)
                else:
                    self._ingroup.append(section)


    def load_files(self):
        for sample in self.samples:
            if isinstance(self.samples[sample].loci, (str, bytes)):
                self.samples[sample].loci = BEDNameMapFile(
                    self.samples[sample].loci
                )
            if isinstance(self.samples[sample].references, (str, bytes)):
                self.samples[sample].references = AssemblyReportFile(
                    self.samples[sample].references,
                    map_assigned_molecule=self.map_assigned_molecule
                )
                
                    
    def from_file(self, infile, **kwargs):
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
            with open(infile, **kwargs) as fp:
                self._parse(fp)

        
    def clear(self):
        self._outgroup = []
        self._ingroup = []
        self.tree = None
        self.samples = dict()
        self.filename = None        

        
    @property
    def species(self):
        return self.samples

    
    @species.setter
    def species(self, samples):
        self.samples = samples


    def has_sample(self, sample_id):
        return sample_id in self.samples


    def has_samples(self, samples):
        return all(map(self.has_sample, samples))
            

    def check_sample(self, sample):
        if not self.has_sample(sample):
            raise KeyError("Sample not found in config: %s" % str(sample))


    def check_samples(self, samples):
         for sample in samples:
             self.check_sample(sample)


             
