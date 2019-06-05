# import utils
from .pySNIDutils import read_output_file, read_lnw

# try to import core functions, but warn if SNID not installed
try:
    from .pySNID import SNID_type, SNID_subtype, SNID_redshift, SNID_age, pySNID
    from .pySNIDutils import exec_SNID
except ImportError as ie:
    print('Warning: {}'.format(ie))
    print('core SNID functions will not be available, but some utility functions are still usable')
