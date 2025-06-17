
'''
Frontend for using logging module as terminal feedback facility.

$Id: feedback.py 3740 2014-10-17 12:39:49Z m221078 $
'''

import logging
import sys

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.WARNING)
logging.addLevelName(logging.DEBUG, 'Debug')
logging.addLevelName(logging.INFO, 'Note')
logging.addLevelName(logging.WARNING, 'Hey')
logging.addLevelName(logging.ERROR, 'Oops')
logging.addLevelName(logging.CRITICAL, 'Sorry')

info = logging.info
warning = logging.warning
error = logging.error
critical = logging.critical

def die(message, *args, **kwargs):
    status = kwargs.pop('status', 1)
    if kwargs:
        (key, dont_care) = kwargs.popitem()
        raise TypeError("die() got an unexpected keyword argument '%s'"%key)
    error(message, *args)
    sys.exit(status)

