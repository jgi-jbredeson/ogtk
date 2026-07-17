#!/usr/bin/env python3

import os
import sys
    
try:
    import og
    import signal
except ImportError as error:
    sys.stderr.write("[%s] %s: %s" % (
        os.path.basename(__file__),
        error.__class__.__name__,
        str(error))
    )
    sys.exit(1)
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

from og import (
    __authors__,
    __contact__,
    __program__,
    __pkgname__,
    __version__
)
__purpose__ = 'Tools for working with orthogroups'

from og.constants import (
    _EMPTY,
)
from og.logging import (
    log_exception,
    warn
)
import og.cli.utils

_MAIN_FLAGS = ('-h', '--help', '-v', '--version')



# def _log_exception(exception, message, traceback):
#     """
# Exception handler hook to redirect errors to be reported 
# by the logging module
#     """
#     logging.error("%s: %s" % (exception.__name__, message))
#     sys.exit(1)
    


def _flag_in(args):
    for flag in _MAIN_FLAGS:
        if flag in args:
            return True
    return False

    
    
def usage(message=None, exitcode=1, stream=sys.stderr):
    message = _EMPTY if message is None else "ERROR: %s\n\n\n" % message

    stream.write("\n")
    stream.write("Program: %s (%s)\n" % (__program__, __purpose__))
    stream.write("Version: %s %s\n" % (__pkgname__, __version__))
    stream.write("Contact: %s\n" % __contact__)
    stream.write("\n")
    stream.write("Usage: %s <command>\n" % __program__)
    stream.write("\n")
    #------------|----+----|----+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #            0        10        20        30        40        50        60        70        80
    stream.write("Commands:\n")
    stream.write("  add-singletons Append missing members to an orthogroups file\n")
    stream.write("  count          Count orthogroups members\n")
    stream.write("  filter         Filter orthogroups files in various ways\n")
    stream.write("  join           Join orthogroups files on member overlap\n")
    stream.write("  score          Use manual clustering to classify orthogroups\n")
    stream.write("  subset         Subset orthogroups file members or group IDs\n")
    stream.write("  update         Reassign groups using posterior probabilities\n")
    stream.write("  upset          Print upset plot describing an orthogroups file\n")
    stream.write("  view           Manipulate/convert orthogroups files\n")
    stream.write("\n")
    stream.write("\n%s" % message)
    
    sys.exit(exitcode)
    
    

if __name__ == '__main__':
    # sys.excepthook = log_exception
    og.cli.utils.MAIN_ENTRY = True
    if len(sys.argv) < 2:
        usage()
    elif _flag_in(sys.argv[1:]):
        usage()
    else:
        main = usage
        command_name = sys.argv[1]
        if command_name == 'add-singletons':
            from og.cli.add_singletons import main
        elif command_name == 'count':
            from og.cli.count import main
        elif command_name == 'filter':
            from og.cli.filter import main
        elif command_name == 'join':
            from og.cli.join import main
        elif command_name == 'score':
            from og.cli.score import main
        elif command_name == 'subset':
            from og.cli.subset import main
        elif command_name == 'update':
            from og.cli.udpate import main
        elif command_name == 'upset':
            from og.cli.upset import main
        elif command_name == 'view':
            from og.cli.view import main
        else:
            usage("Unrecognized command: %s\n" % command_name)

        from warnings import catch_warnings
        with catch_warnings(record=True) as stack:
            main(sys.argv[2:])
            for warning in stack:
                warn(str(warning.message))
                

