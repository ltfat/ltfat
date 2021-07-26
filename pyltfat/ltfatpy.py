# -*- coding: utf-8 -*-
# This file is an extension of oct2py (Copyright (c) oct2py developers, distributed under the MIT License),
# to simplify calling the LTFAT from Python.
# Distributed under the terms of the GPL v3 License.

import logging
import os.path as osp

import numpy as np
from metakernel.pexpect import EOF, TIMEOUT

from oct2py.core import Oct2Py
from oct2py.io import read_file, write_file, Cell
from oct2py.utils import Oct2PyError
from oct2py.compat import string_types
from oct2py.dynamic import OctavePtr



class LTFAT(Oct2Py):

    """Manages an Octave session.

    Uses MAT files to pass data between Octave and Numpy.
    The function must either exist as an m-file in this directory or
    on Octave's path.
    The first command will take about 0.5s for Octave to load up.
    The subsequent commands will be much faster.

    You may provide a logger object for logging events, or the oct2py.get_log()
    default will be used.  When calling commands, logger.info() will be used
    to stream output, unless a `stream_handler` is provided.

    Parameters
    ----------
    logger : logging object, optional
        Optional logger to use for Oct2Py session
    timeout : float, optional
        Timeout in seconds for commands
    oned_as : {'row', 'column'}, optional
        If 'column', write 1-D numpy arrays as column vectors.
        If 'row', write 1-D numpy arrays as row vectors.}
    temp_dir : str, optional
        If specified, the session's MAT files will be created in the
        directory, otherwise a default directory is used.  This can be
        a shared memory (tmpfs) path.
    convert_to_float : bool, optional
        If true, convert integer types to float when passing to Octave.
    backend: string, optional
        The graphics_toolkit to use for plotting.
    """

    def __init__(self, logger=None, timeout=None,
                 oned_as='row', temp_dir=None, convert_to_float=True,
                 backend=None):
        """Start Octave and set up the session.
        """
        Oct2Py.__init__(self)


    def _feval(self, func_name, func_args=(), dname='', nout=0,
              timeout=None, stream_handler=None, store_as='', plot_dir=None):
        """Run the given function with the given args.
        """
        engine = self._engine
        if engine is None:
            raise Oct2PyError('Session is closed')

        # Set up our mat file paths.
        out_file = osp.join(self.temp_dir, 'writer.mat')
        out_file = out_file.replace(osp.sep, '/')
        in_file = osp.join(self.temp_dir, 'reader.mat')
        in_file = in_file.replace(osp.sep, '/')

        func_args = list(func_args)
        ref_indices = []
        for (i, value) in enumerate(func_args):
            if isinstance(value, OctavePtr):
                ref_indices.append(i + 1)
                func_args[i] = value.address
        ref_indices = np.array(ref_indices)

        # Save the request data to the output file.
        req = dict(func_name=func_name, func_args=tuple(func_args),
                   dname=dname or '', nout=nout,
                   store_as=store_as or '',
                   ref_indices=ref_indices)

        write_file(req, out_file, oned_as=self._oned_as,
                   convert_to_float=self.convert_to_float)

        # Set up the engine and evaluate the `_pyeval()` function.
        engine.line_handler = stream_handler or self.logger.info
        if timeout is None:
            timeout = self.timeout

        try:
            engine.eval('_pyeval("%s", "%s");' % (out_file, in_file),
                        timeout=timeout)
        except KeyboardInterrupt as e:
            stream_handler(engine.repl.interrupt())
            raise
        except TIMEOUT:
            stream_handler(engine.repl.interrupt())
            raise Oct2PyError('Timed out, interrupting')
        except EOF:
            if not self._engine:
                return
            stream_handler(engine.repl.child.before)
            self.restart()
            raise Oct2PyError('Session died, restarting')

        # Read in the output.
        resp = read_file(in_file, self)
        print(resp)
        
        if resp['err']:
            msg = self._parse_error(resp['err'])
            raise Oct2PyError(msg)

        result = resp['result'].ravel().tolist()
        if isinstance(result, list) and len(result) == 1:
            result = result[0]

        # Check for sentinel value.
        if (isinstance(result, Cell) and
                result.size == 1 and
                isinstance(result[0], string_types) and
                result[0] == '__no_value__'):
            result = None

        if plot_dir:
            self._engine.make_figures(plot_dir)

        return result
