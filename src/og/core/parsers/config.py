
import yaml
from math import inf as _POS_INF
from og.core.parsers.bed import BEDNameMap
from og.core.parsers.newick import IntervalNewickTree
from og.core.parsers.assembly_report import AssemblyReport
from og.constants import _PYTHON_VERSION
from og.core.io import open, is_stream

if _PYTHON_VERSION < (3,7):
    from collections import OrderedDict as dict

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


class _SpeciesRecord(object):
    def __init__(self, species=None, loci=None, references=None, unplaced_id=None, is_outgroup=False):
        self.species = species
        self.loci = loci
        self.references = references
        self.unplaced_id = unplaced_id
        self.is_outgroup = is_outgroup
        

class SpeciesConfig(object):
    def __init__(self, infile, load_files=False):
        self.clear()
        self.filename = None
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
            else:  # a species ID
                self.species[section] = _SpeciesRecord(species=section)
                if 'loci' in yamldata[section]:
                    self.species[section].loci = yamldata[section]['loci']  # BEDNameMap(yamldata[section]['loci'])
                else:
                    raise ValueError('`loci` key not defined for %s' % section)
                
                if 'references' in yamldata[section]:
                    self.species[section].references = yamldata[section]['references'] # AssemblyReport(yamldata[section]['references'])
                else:
                    raise ValueError('`references` key not defined for %s' % section)
                
                if 'unplaced_id' in yamldata[section]:
                    self.species[section].unplaced_id = yamldata[section]['unplaced_id']
                else:
                    raise ValueError('`unplaced_id` key not defined for %s' % section)
                
                if 'outgroup' in yamldata[section]:
                    self.species[section].is_outgroup = to_bool(str(yamldata[section]['outgroup']))

                if self.species[section].is_outgroup:
                    self._outgroup.append(section)
                else:
                    self._ingroup.append(section)


    def load_files(self):
        for species in self.species:
            if isinstance(self.species[species].loci, (str, bytes)):
                self.species[species].loci = BEDNameMap(self.species[species].loci)
            if isinstance(self.species[species].references, (str, bytes)):
                self.species[species].references = AssemblyReport(self.species[species].references)
                
                    
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
        self.species = dict()
        self.filename = None        
