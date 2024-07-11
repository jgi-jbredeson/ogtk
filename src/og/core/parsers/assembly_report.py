
from og.constants import (
    _PYTHON_VERSION,
    _COMMENT,
    _TAB
)
from ..io import open, is_stream


if _PYTHON_VERSION < (3,7):
    import collections.OrderedDict as dict
    

num = len
ROLE_ASM_MOL = 0x1
ROLE_UNLOC_SCAF = 0x2
ROLE_UNPLC_SCAF = 0x4
ROLE_UNLOC_CTG = 0x8
ROLE_UNPLC_CTG = 0x10
ROLE_ALT_SCAF = 0x20
ROLE_ALT_CTG = 0x40
ROLE_FIX_PATCH = 0x80
ROLE_NOVEL_PATCH = 0x100
ROLE_IGNORED = 0x200

UNIT_PRIMARY = 0x1
UNIT_ALTERNATE = 0x2
UNIT_PATCHES = 0x4
UNIT_NONNUCLEAR = 0x8

TYPE_CHROMOSOME = 0x1


class AssemblyReportFormatError(Exception):
    pass


class _AssemblyReportRecord(object):
    def __init__(self,
                 sequence_name, sequence_length=-1, sequence_role=0,
                 assigned_molecule=-1, assigned_type=0, genbank_accession=None,
                 refseq_accession=None, assembly_unit=0, ucsc_name=None):
        self.sequence_name = sequence_name
        self.sequence_role = sequence_role
        self.sequence_length = sequence_length
        self.assigned_molecule = assigned_molecule
        self.assigned_type = assigned_type
        self.genbank_accession = genbank_accession
        self.refseq_accession = refseq_accession
        self.assembly_unit = assembly_unit
        self.ucsc_name = ucsc_name
        self.is_chr = False
        
    
class AssemblyReport(dict):
    def __init__(self, infile=None, **kwargs):
        dict.__init__(self)
        if infile is not None:
            self.from_file(infile, **kwargs)
            
        
    def _read_file(self, infile):
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
                record.sequence_role |= ROLE_ALT_SCAF
                record.assembly_unit |= UNIT_ALTERNATE
            elif sequence_role == 'alt-contig':
                record.sequence_role |= ROLE_ALT_CTG
                record.assembly_unit |= UNIT_ALTERNATE
            elif sequence_role == 'fix-patch':
                record.sequence_role |= ROLE_FIX_PATCH
            elif sequence_role == 'novel-patch':
                record.sequence_role |= ROLE_NOVEL_PATCH
            elif sequence_role == 'assembled-molecule':
                record.sequence_role |= ROLE_ASM_MOL
            elif sequence_role == 'unlocalized-scaffold':
                record.sequence_role |= ROLE_UNLOC_SCAF
            elif sequence_role == 'unplaced-scaffold':
                record.sequence_role |= ROLE_UNPLC_SCAF
            elif sequence_role == 'unlocalized-contig':
                record.sequence_role |= ROLE_UNLOC_CTG
            elif sequence_role == 'unplaced-contig':
                record.sequence_role |= ROLE_UNPLC_CTG

            if assigned_type == 'chromosome':
                record.assigned_type |= TYPE_CHROMOSOME
                record.assigned_molecule = fields[2]
                
            if assembly_unit == 'primary assembly':
                record.assembly_unit |= UNIT_PRIMARY
            elif assembly_unit == 'non-nuclear':
                record.assembly_unit |= UNIT_NONNUCLEAR
            elif assembly_unit == 'patches':
                record.assembly_unit |= UNIT_PATCHES
                
            if fields[4].lower() != 'na':
                record.genbank_accession = fields[4]
            if fields[6].lower() != 'na':
                record.refseq_accession = fields[6]
            if fields[9].lower() != 'na':
                record.ucsc_name = fields[9]

            if ((record.assembly_unit & UNIT_PRIMARY) and
                (record.sequence_role & ROLE_ASM_MOL) and
                (record.assigned_type & TYPE_CHROMOSOME)):
                record.is_chr = True
                
            self[record.sequence_name] = record

    
    def from_file(self, infile, **kwargs):
        self.clear()
        if is_stream(infile):
            self._read_file(infile)
        else:
            if 'mode' in kwargs:
                if 'a' in kwargs['mode'] or \
                   'w' in kwargs['mode']:
                    raise ValueError("%s() constructor is read-only" % (
                        self.__class__.__name__
                    ))
            else:
                kwargs['mode'] = 'rt'
            with open(infile, **kwargs) as fd:
                self._read_file(fd)


    def from_string(self, instring):
        import io
        self.clear()
        self._read_file(io.StringIO(instring))

        
