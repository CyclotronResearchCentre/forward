###
# Sample data loading code. Adapted from MNE Python.
#
# All credit to:
##
####
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#          Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#          Eric Larson <larson.eric.d@gmail.com>
# License: BSD Style.

from .utils import _data_path

def data_path(path=None, force_update=False, update_path=True,
              name="example", download=True, verbose=None):
    return _data_path(path=path, force_update=force_update,
                      update_path=update_path, name=name,
                      download=download,
                      verbose=verbose)