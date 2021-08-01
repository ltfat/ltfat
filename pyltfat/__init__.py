# -*- coding: utf-8 -*-
# This file is an extension of oct2py (Copyright (c) oct2py developers, distributed under the MIT License),
# to simplify calling the LTFAT from Python.
# Distributed under the terms of the GPL v3 License.

from .ltfatpy import LTFAT
from oct2py.utils import get_log, Oct2PyError
import os

dirname = os.getcwd()
foldername = os.path.join(os.sep, dirname, 'ltfat', 'pyltfat', 'files')
#os.chdir(foldername)

__all__ = ['ltfatpy']

#import logging
#logging.basicConfig(level=logging.DEBUG)
try:
    ltfatpy = LTFAT(temp_dir=foldername)
except Oct2PyError as e:
    print(e)

#ltfatpy.run('ltfatstart')
#os.chdir(dirname)

def kill_octave():
    """Kill all octave instances (cross-platform)."""
    import os
    if os.name == 'nt':
        os.system('taskkill /im octave /f')
    else:
        os.system('killall -9 octave')
        os.system('killall -9 octave-cli')
    ltfatpy.restart()
