
import sys
import signal
import logging

# from logging import info, debug
# from logging import warning as warn

from og import __program__
from og.constants import _SPACE

logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="[%(name)s] %(levelname)s: %(message)s"
)
signal.signal(signal.SIGPIPE,signal.SIG_DFL)


logger = logging.getLogger(__program__)



def log_exception(exception, message, traceback):
    """
Exception handler hook to redirect errors to be reported 
by the logging module
    """
    logging.error("%s: %s" % (exception.__name__, message))
    sys.exit(1)

    
    
def set_level(level):
    logger.setLevel(level)



def info(*args):
    logger.info(_SPACE.join(args))


    
def warn(*args):
    logger.warning(_SPACE.join(args))

    

def debug(*args):
    logger.debug(_SPACE.join(args))

    
