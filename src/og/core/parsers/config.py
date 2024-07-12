
from configparser import ConfigParser
from og.core.parsers.bed import BEDNameMap
from og.core.parsers.newick import IntervalNewickTree
from og.core.parsers.assembly_report import AssemblyReport
from og.constants import (
    _PYTHON_VERSION,
)
from ..io import open, is_stream

if _PYTHON_VERSION < (3,7):
    import collections.OrderedDict as dict

        
class SpeciesConfig(object):
    def __init__(self, infile):
        self.clear()
        self.filename = None
        if infile is not None:
            self.from_file(infile)

            
    def _parse(self, infile):
        parser = ConfigParser()
        parser.read_file(infile)
        for section in parser.sections():
            if section == 'tree':
                self.tree = dict(parser[section])
                self.tree['ploidy'] = IntervalNewickTree(self.tree['ploidy'])
                
            else:  # a species ID
                self.species[section] = {}
                if 'loci' in parser[section]:
                    self.species[section]['loci'] = BEDNameMap(parser[section]['loci'])
                else:
                    raise ValueError('`loci` key not provided for %s' % section)
                if 'references' in parser[section]:
                    self.species[section]['references'] = AssemblyReport(parser[section]['references'])
                else:
                    raise ValueError('`references` key not provided for %s' % section)
                
                if 'outgroup' in parser[section] and parser[section].getboolean('outgroup'):
                    self.species[section]['outgroup'] = True
                    self._outgroup.append(section)
                else:
                    self.species[section]['outgroup'] = False
                    self._ingroup.append(section)

                    
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
