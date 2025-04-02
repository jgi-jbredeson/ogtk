
import sys

__version__ = '__PACKAGE_VERSION__'
__authors__ = '__PACKAGE_AUTHORS__'
__contact__ = '__PACKAGE_CONTACT__'

__all__ = (
    '_PYTHON_VERSION',
    '_EXCLAMATION_MARK', '_XMARK', '_BANG',
    '_BACK_SLASH','_BSLASH',
    '_FORWARD_SLASH', '_FWRDSLASH', '_FSLASH',
    '_COLON',
    '_COMMA',
    '_COMMENT', '_HASH',
    '_DASH', '_MINUS',
    '_PLUS',
    '_EMPTY',
    '_EOL',
    '_EQUAL',
    '_PERIOD','_DOT',
    '_LEFT_PARENTHESIS',  '_LPAREN',
    '_RIGHT_PARENTHESIS', '_RPAREN',
    '_SEMICOLON',
    '_DOUBLEQUOTE', '_2QUOTE',
    '_SINGLEQUOTE', '_1QUOTE',
    '_QUESTION_MARK', '_QMARK',
    '_PIPE',
    '_SPACE',
    '_TAB',
    'range',
    'dict'
)


_PYTHON_VERSION = sys.version_info[:2]

_EXCLAMATION_MARK = _XMARK = _BANG = '!'
_BACK_SLASH = _BSLASH = '\\'
_FORWARD_SLASH = _FWRDSLASH = _FSLASH = '/'
_COLON = ':'
_COMMA = ','
_COMMENT = _HASH = '#'
_DASH = _MINUS = '-'
_PLUS = '+'
_EMPTY = ''
_EOL = '\n'
_EQUAL = '='
_PERIOD = _DOT = '.'
_LEFT_PARENTHESIS = _LPAREN = '('
_RIGHT_PARENTHESIS = _RPAREN= ')'
_SEMICOLON = ';'
_DOUBLEQUOTE = _2QUOTE = '"'
_SINGLEQUOTE = _1QUOTE = "'"
_QUESTION_MARK = _QMARK = '?'
_PIPE = '|'
_SPACE = ' '
_TAB = '\t'

_dummy = range
if _PYTHON_VERSION < (3,):
    range = xrange
else:
    range = _dummy

_dummy = dict
if _PYTHON_VERSION < (3,7):
    from collections import OrderedDict as dict
else:
    dict = _dummy
