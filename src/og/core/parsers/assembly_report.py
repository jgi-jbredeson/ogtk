
import sys
from og.constants import (
    _COMMENT,
    _SPACE,
    _EOL,
    _TAB,
    dict
)
from og.core.io import open, is_stream


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

REL_EQUAL = 0x2000
REL_CMP = 0x4000

_VALID_COLUMNS = (
    'Sequence-Name',
    'Sequence-Role',
    'Assigned-Molecule',
    'Assigned-Molecule-Location/Type',
    'GenBank-Accn',
    'Relationship',
    'RefSeq-Accn',
    'Assembly-Unit',
    'Sequence-Length',
    'UCSC-style-name'
)
_VALID_ATTRIBUTES = (
    'sequence_name',
    'sequence_role',
    'assigned_molecule',
    'assigned_type',
    'genbank_accession',
    'relationship',
    'refseq_accession',
    'assembly_unit',
    'sequence_length',
    'ucsc_name'
)

def NoneToNa(value):
    return 'na' if value is None else value


class AssemblyReportFormatError(Exception):
    pass


class _AssemblyReportRecord(object):
    def __init__(self,
                 sequence_name, sequence_length=-1,
                 assigned_molecule=None, assigned_type=None,
                 assembly_unit=None, flags=0,
                 genbank_accession=None, refseq_accession=None, ucsc_name=None):
        self.sequence_name = sequence_name
        self.sequence_length = sequence_length
        self.flags = flags
        self.assigned_molecule = assigned_molecule
        self.assigned_type = assigned_type
        self.assembly_unit = assembly_unit
        self.genbank_accession = genbank_accession
        self.refseq_accession = refseq_accession
        self.ucsc_name = ucsc_name

    def __repr__(self):
        return '%s(%s)' % (
            self.__class__.__name__,
            ', '.join([
                '%s=%s' % (
                    attr,
                    str(getattr(self, attr, None))) \
                    for attr in _VALID_ATTRIBUTES
            ])
        )

    def __str__(self):
        return _TAB.join((
            str(NoneToNa(getattr(self, attr, None))) for attr in _VALID_ATTRIBUTES
        ))
    
    @property
    def sequence_role(self):
        if self.flags & ROLE_ASM_MOL:
            return 'assembled-molecule'
        elif self.flags & ROLE_UNPLC_SCAF:
            return 'unplaced-scaffold'
        elif self.flags & ROLE_UNPLC_CTG:
            return 'unplaced-contig'
        elif self.flags & ROLE_UNLOC_SCAF:
            return 'unlocalized-scaffold'
        elif self.flags & ROLE_UNLOC_CTG:
            return 'unlocalized-contig'
        elif self.flags & ROLE_ALT_SCAF:
            return 'alt-scaffold'
        elif self.flags & ROLE_ALT_CTG:
            return 'alt-contig'
        elif self.flags & ROLE_FIX_PATCH:
            return 'fix-patch'
        elif self.flags & ROLE_NOVEL_PATCH:
            return 'novel-patch'
        return None

    @property
    def relationship(self):
        if self.flags & REL_EQUAL:
            return '='
        elif self.flags & REL_CMP:
            return '<>'
        else:
            return None
    
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

    
    
class AssemblyReportFile(dict):
    def __init__(self, infile=None, **kwargs):
        dict.__init__(self)
        self.filename = None
        self.map_assigned_molecule = False
        if 'map_assigned_molecule' in kwargs:
            self.map_assigned_molecule = kwargs['map_assigned_molecule']
            del(kwargs['map_assigned_molecule'])
        if infile is not None:
            self.from_file(infile, **kwargs)
            
        
    def _parse(self, infile):
        line_count = 0
        assembled_molecules = dict()
        for line in infile:
            line = line.lstrip().rstrip('\r\n')
            line_count += 1
            
            if line.startswith(_COMMENT):
                continue

            fields = line.split(_TAB)

            if num(fields) < 10:
                raise AssemblyReportFormatError(
                    "Ten fields expected, line %d" % line_count
                )

            for i in range(num(fields)):
                if fields[i].lower() == "na":
                    fields[i] = None
        
            try:
                fields[8] = int(fields[8])
            except ValueError:
                raise AssemblyReportFormatError(
                    "Numeric field expected, line %d column %d" % (
                        line_count, 9
                    )) from None
            
            record = _AssemblyReportRecord(
                sequence_name=fields[0],
                assigned_molecule=fields[2],
                assigned_type=fields[3],
                genbank_accession=fields[4],
                refseq_accession=fields[6],
                assembly_unit=fields[7],
                sequence_length=fields[8],
                ucsc_name=fields[9]
            )

            sequence_role = None \
                if   fields[1] is None \
                else fields[1].lower()
            assigned_type = None \
                if   record.assigned_type is None \
                else record.assigned_type.lower()
            assembly_unit = None \
                if   record.assembly_unit is None \
                else record.assembly_unit.lower()
            relationship = fields[5]
            
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
                assembled_molecules[record.assigned_molecule] \
                    = record.sequence_name
                record.flags |= ROLE_ASM_MOL
            elif sequence_role == 'unlocalized-scaffold':
                record.flags |= ROLE_UNLOC_SCAF
            elif sequence_role == 'unplaced-scaffold':
                record.flags |= ROLE_UNPLC_SCAF
            elif sequence_role == 'unlocalized-contig':
                record.flags |= ROLE_UNLOC_CTG
            elif sequence_role == 'unplaced-contig':
                record.flags |= ROLE_UNPLC_CTG

            # if assigned_type == 'chromosome':
            #     record.assigned_molecule = fields[2]
                
            if assembly_unit == 'primary assembly':
                record.flags |= UNIT_PRIMARY
            elif assembly_unit == 'non-nuclear':
                record.flags |= UNIT_NONNUCLEAR
            elif assembly_unit == 'patches':
                record.flags |= UNIT_PATCHES

            if relationship == '=':
                record.flags |= REL_EQUAL
            elif relationship == '<>':
                record.flags |= REL_CMP
                
            # if fields[4] is not None:
            #     record.genbank_accession = fields[4]
            # if fields[6] is not None:
            #     record.refseq_accession = fields[6]
            # if fields[9] is not None:
            #     record.ucsc_name = fields[9]

            self[record.sequence_name] = record

        # Map Sequence-Name of assembled-molecule onto Assigned-Molecule field:
        if self.map_assigned_molecule:
            for record in self.values():
                if record.assigned_molecule is None:
                    continue
                elif record.assigned_molecule in assembled_molecules:
                    record.assigned_molecule \
                        = assembled_molecules[record.assigned_molecule]
                else:
                    raise KeyError("Assigned molecule: %s" % record.assigned_molecule)


    def _format_header(self):
        return _COMMENT + _SPACE + _TAB.join(_VALID_COLUMNS)

    
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

        stream.write(self._format_header() + _EOL)
        for record in self.values():
            stream.write(str(record) + _EOL)

        if close:
            stream.close()

                
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


