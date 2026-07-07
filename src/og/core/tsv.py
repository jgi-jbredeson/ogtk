
from og.constants import (
    _COMMENT,
    _EMPTY,
    _EOL,
    _TAB,
    dict
)
from og.core.compressio import open

num = len

def read_to_dict(filename, key=0, fields=(1,), sep=None):
    data_dict = dict()
    max_fields = max(fields)
    num_lines = 0
    with open(filename,'r') as tsv_file:
        for line in tsv_file:
            line = line.rstrip(_EOL)

            num_lines += 1
            
            if line.strip() == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            _fields = line.split(sep)

            if num(_fields) < max_fields:
                raise IndexError(
                    "Too few fields, at least %d expected, line %d" % (
                        max_fields, num_lines
                    )
                )

            if _fields[key] in data_dict:
                raise KeyError("Duplicate key: %s" % _fields[key])
            else:
                data_dict[_fields[key]] = []
                for field in fields:
                    data_dict[_fields[key]].append(_fields[field])

    return data_dict


def read_to_list(filename, field=0, sep=None):
    data_list = list()
    max_fields = field + 1
    num_lines = 0
    with open(filename,'r') as tsv_file:
        for line in tsv_file:
            line = line.rstrip(_EOL)

            num_lines += 1

            if line.strip() == _EMPTY or \
               line.startswith(_COMMENT):
                continue

            _fields = line.split(sep)

            if num(_fields) < max_fields:
                raise IndexError(
                    "Too few fields, at least %d expected, line %d" % (
                        max_fields, num_lines
                    )
                )

            data_list.append(_fields[field])
            
    return data_list

