
from og.constants import (
    _PYTHON_VERSION,
    _COMMENT,
    _TAB
)
from og.core.io import open, is_stream


if _PYTHON_VERSION < (3,7):
    import collections.OrderedDict as dict
    

num = len

ROLE_ASM_MOL = 0x1
ROLE_UNPLC_SCAF = 0x2
ROLE_UNPLC_CTG = 0x4
ROLE_UNLOC_SCAF = 0x8
ROLE_UNLOC_CTG = 0x10
ROLE_ALT_SCAF = 0x20
ROLE_ALT_CTG = 0x40
ROLE_FIX_PATCH = 0x80
ROLE_NOVEL_PATCH = 0x100

UNIT_PRIMARY = 0x200
UNIT_ALTERNATE = 0x400
UNIT_PATCHES = 0x800
UNIT_NONNUCLEAR = 0x1000


class AssemblyReportFormatError(Exception):
    pass


class _AssemblyReportRecord(object):
    def __init__(self,
                 sequence_name, sequence_length=-1, flags=0,
                 assigned_molecule=None, assigned_type=None,
                 genbank_accession=None, refseq_accession=None, ucsc_name=None):
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
        self.flags = flags
        self.assigned_molecule = assigned_molecule
        self.assigned_type = assigned_type        
        self.genbank_accession = genbank_accession
        self.refseq_accession = refseq_accession
        self.ucsc_name = ucsc_name
        
    @property
    def is_primary(self):
        return bool(self.flags & UNIT_PRIMARY)

    @property
    def is_alternate(self):
        return bool(self.flags & (
            UNIT_ALTERNATE | ROLE_ALT_SCAF | ROLE_ALT_CTG
        ))

    @property
    def is_patch(self):
        return bool(self.flags & (
            UNIT_PATCHES | ROLE_FIX_PATCH | ROLE_NOVEL_PATCH
        ))

    @property
    def is_nonnuclear(self):
        return bool(self.flags & UNIT_NONNUCLEAR)

    @property
    def is_assembled_molecule(self):
        return bool(self.flags & ROLE_ASM_MOL)

    @property
    def is_unplaced(self):
        return bool(self.flags & (ROLE_UNPLC_CTG | ROLE_UNPLC_SCAF))

    @property
    def is_unlocalized(self):
        return bool(self.flags & (ROLE_UNLOC_CTG | ROLE_UNLOC_SCAF))
    
    @property
    def is_placed(self):
        return not self.is_unplaced

    @property
    def is_localized(self):
        return not (self.is_unplaced or self.is_unlocalized)

    
    
class AssemblyReport(dict):
    def __init__(self, infile=None, **kwargs):
        dict.__init__(self)
        self.filename = None
        if infile is not None:
            self.from_file(infile, **kwargs)
            
        
    def _parse(self, infile):
        line_count = 0
        for line in infile:
            line = line.lstrip().rstrip('\r\n')
            line_count += 1
            
            if line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)

            if num(fields) != 10:
                raise AssemblyReportFormatError(
                    "Ten fields expected, line %d" % line_count
                )

            try:
                fields[8] = int(fields[8])
            except ValueError:
                raise AssemblyReportFormatError(
                    "Numeric field expected, line %d column %d" % (
                        line_count, 9
                    )) from None
                
            record = _AssemblyReportRecord(fields[0], fields[8])

            sequence_role = fields[1].lower()
            assigned_type = fields[3].lower()
            assembly_unit = fields[7].lower()
            if sequence_role == 'alt-scaffold':
                record.flags |= ROLE_ALT_SCAF
                record.flags |= UNIT_ALTERNATE
            elif sequence_role == 'alt-contig':
                record.flags |= ROLE_ALT_CTG
                record.flags |= UNIT_ALTERNATE
            elif sequence_role == 'fix-patch':
                record.flags |= ROLE_FIX_PATCH
            elif sequence_role == 'novel-patch':
                record.flags |= ROLE_NOVEL_PATCH
            elif sequence_role == 'assembled-molecule':
                record.flags |= ROLE_ASM_MOL
            elif sequence_role == 'unlocalized-scaffold':
                record.flags |= ROLE_UNLOC_SCAF
            elif sequence_role == 'unplaced-scaffold':
                record.flags |= ROLE_UNPLC_SCAF
            elif sequence_role == 'unlocalized-contig':
                record.flags |= ROLE_UNLOC_CTG
            elif sequence_role == 'unplaced-contig':
                record.flags |= ROLE_UNPLC_CTG

            if assigned_type == 'chromosome':
                record.assigned_molecule = fields[2]
                
            if assembly_unit == 'primary assembly':
                record.flags |= UNIT_PRIMARY
            elif assembly_unit == 'non-nuclear':
                record.flags |= UNIT_NONNUCLEAR
            elif assembly_unit == 'patches':
                record.flags |= UNIT_PATCHES
                
            if fields[4].lower() != 'na':
                record.genbank_accession = fields[4]
            if fields[6].lower() != 'na':
                record.refseq_accession = fields[6]
            if fields[9].lower() != 'na':
                record.ucsc_name = fields[9]

            self[record.sequence_name] = record

    
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


    def from_string(self, instring):
        import io
        self.clear()
        self._parse(io.StringIO(instring))

        
    def clear(self):
        dict.clear(self)
        self.filename = None



def is_placed(record):
    return not (record.is_unplaced or record.assigned_molecule is None)


def is_chr(record):
    return is_placed(record) and record.is_assembled_molecule


