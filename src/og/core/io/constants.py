
import sys

_PYTHON_VERSION = sys.version_info[:2]

STDIO = '-'

COMPRESSION_SUFFIX = {
    'lzma'  : ('.lzma','.xz'),
    'gzip'  : ('.gzip','-gz','.z','-z','_z'),
    'bzip2' : ('.bzip2','.bz2','.bz'),
    'bgzip' : ('.bgzip','.bgzf','.bgz','.gz')
}

COMPRESSION_NAME_MAP = {
    'lzma'    : 'lzma',
    'xz'      : 'lzma',
    'gzip'    : 'gzip',
    'gz'      : 'gzip',
    'bzip2'   : 'bzip2',
    'bz2file' : 'bzip2',
    'bz2'     : 'bzip2',
    'bgzip'   : 'bgzip',
    'bgz'     : 'bgzip'
}
