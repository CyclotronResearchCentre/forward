###
# Sample data loading code. Adapted from MNE Python.
#
# All credit to:
##
####
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#          Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#          Eric Larson <larson.eric.d@gmail.com>
#          Denis Egnemann <d.engemann@fz-juelich.de>
# License: BSD Style.

import os
import os.path as op
import shutil
import tarfile
import logging
logger = logging.getLogger('forwardlog')

from mne.utils import _fetch_file

def _data_path(path=None, force_update=False, update_path=True,
               download=True, name=None, verbose=None):
    if path is None:
        path = op.join(os.environ["FWD_DIR"],"examples")

    if not isinstance(path, basestring):
        raise ValueError('path must be a string or None')

    if name.lower() == 'example':
        archive_name = "ForwardSample.tar.gz"
        url = "https://www.dropbox.com/s/dr7qx6hbv7myxdc/" + archive_name + "?dl=1"
        folder_name = "ForwardSample"
        folder_path = op.join(path, folder_name)
        rm_archive = False
    elif name.lower() == 'leadfield':
        archive_name = "LeadfieldSample.tar.gz"
        url = "https://www.dropbox.com/s/3qic6ma3umzp2jg/" + archive_name + "?dl=1"
        folder_name = "LeadfieldSample"
        folder_path = op.join(path, folder_name)
        rm_archive = False
    elif name.lower() == 'simnibs':
        archive_name = "simnibs_example.tar.gz"
        url = "http://simnibs.org/_media/" + archive_name
        folder_name = "simnibs_example"
        folder_path = op.join(path, folder_name)
        rm_archive = False
    else:
        raise ValueError('Sorry, the dataset "%s" is not known.' % name)

    if not op.exists(folder_path) and not download:
        return ''

    if not op.exists(folder_path) or force_update:
        logger.info('Sample data archive %s not found at:\n%s\n'
                    'It will be downloaded and extracted at this location.'
                    % (archive_name, folder_path))

        archive_name = op.join(path, archive_name)
        rm_archive = True
        if op.exists(archive_name):
            msg = ('Archive already exists at %r. Overwrite it '
                   '(y/[n])? ' % archive_name)
            answer = raw_input(msg)
            if answer.lower() == 'y':
                os.remove(archive_name)
            else:
                raise IOError('Archive file already exists at target '
                              'location %r.' % archive_name)

        _fetch_file(url, archive_name, print_destination=False)

        if op.exists(folder_path):
            shutil.rmtree(folder_path)

        logger.info('Decompressing the archive: ' + archive_name)
        logger.info('... please be patient, this can take some time')
        for ext in ['gz', 'bz2']:  # informed guess (and the only 2 options).
            try:
                tarfile.open(archive_name, 'r:%s' % ext).extractall(path=path)
            except tarfile.ReadError, err:
                logger.info('%s is %s trying "bz2"' % (archive_name, err))

        if rm_archive:
            os.remove(archive_name)

    path = op.abspath(path)
    path = op.join(path, folder_name)

    return path
