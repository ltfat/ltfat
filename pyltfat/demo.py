# This file is an extension of oct2py (Copyright (c) oct2py developers, distributed under the MIT License),
# to simplify calling the LTFAT from Python.
# Distributed under the terms of the GPL v3 License.

# import the LTFAT class
from ltfatpy import LTFAT

# instantiate the lp object that bridges to Octave
lp = LTFAT()
# in Octave, add LTFATs topfolder
lp.addpath('../../ltfat')
# run (the script) ltfatstart.m to set up all paths within Octave
lp.run('ltfatstart')
# optionally, if you want to use the C-backend too
# lp.run('ltfatmex')


# push a list y to the octave workspace
# Integer type arguments will be converted to floating point unless `convert_to_float=False`.
y = [1, 2]
lp.push('y', y)

# get y back into python
lp.pull('y')
#expected: array([[1., 2.]])
# ...you can also get is as a pointer...
ptr = lp.get_pointer('y')
ptr.value
ptr.address
# ...finally, you can also pass it back as an argument to the Octave workspace
# (disp is a built-in Octave function for displaying values of variables)
lp.disp(ptr)

#access the help for LTFAT functions (e.g. the function frame.m)
help(lp.frame)

#run an (ltfat-)function in Octave and return the result to Python
x = lp.feval('svd', lp.hilb(3))
x
# specify three return values
(u, v, d) = lp.feval('svd', lp.hilb(3), nout=3)
u.shape

# Evaluate an Octave command or commands.
lp.eval('for i = 1:3; disp(i);end', stream_handler=lines.append)

# the whole procedure for the frame object
lp.eval("F=frame('wmdct','gauss',40);")
F = lp.frame('wmdct','gauss',40);
lp.feval('frame', 'wmdct', 'gauss', 40)
ptr = lp.get_pointer('F')


