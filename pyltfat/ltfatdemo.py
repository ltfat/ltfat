import os
from ltfat.pyltfat import ltfatpy as lp
# 1. THE SETUP (you need to have the ltfat repository in a subfolder to this script)__________
#tell Octave to add the ltfat path
lp.addpath('./ltfat')
# this could be done in Octave, via lp.addpath('./ltfat') too, but depending on your setup,
# this may give ugly warnings when running ltfatstart - I try to make it less verbose
lp.eval('ltfatstart')

# 2. BASIC FUNCTIONALITY_____________________________________________________________________
#pull is really just for looking at a variable
#push is for specifying variables in python and adding them to the octave workspace

#examples for evaluating ltfat functions:
# 1. evaluation in the Octave workspace:
lp.eval("F=frame('wmdct','gauss',40);")
#F is evaluated in the octave workspace, to access it in Python as a struct, type
lp.pull('F') #this will be verbose
x = lp.pull('F') #this will not be verbose
#...but you can also access it as a pointer...
ptr = lp.get_pointer('F') #this will not be verbose
ptr.value #this will be verbose
ptr.address #this will be verbose

# Specific remarks regarding frames:

# in general, frames are passed as an oct2py.io.Struct
type(x)
# you can access their methods...
x.lengthcoef
# ...but they are stored as strings...
type(x.lengthcoef)
# ...because there is no simple way to parse them to a Python-lambda function directly
# you would have to do this manually, e.g.
x.lengthcoef = lambda Ncoef : Ncoef/lp.framered(x)
x.lengthcoef(100)
# next steps could be automating this

# 2. the elegant version: evaluation in Python (=execute the ltfat function in Octave and return the result to Python)
lp.frame('wmdct','gauss',40) #this will be verbose
F = lp.frame('wmdct','gauss',40) #this will not be verbose
# again, F will be of type oct2py.io.Struct, filled with strings

# 3. evaluation as a function in the Octave workspace (the same as in 1.) applies, but the calling syntax differs)
lp.feval('frame', 'wmdct', 'gauss', 40)


# 3. LTFAT specific Demos:___________________________________________________________________________________
# 1.) plotting a dgt (very elegant: nearly everything happens in Python)
a=10
M=40
L=a*M
h=lp.pherm(L,4)# 4th order hermite function.
c=lp.dgt(h,'gauss',a,M); #semicolon or not makes no difference in Python

# Specific remarks on plotting:

# a) plot using matplotlib in Python (preferred, the syntax is the same as in Octave - Octave internally also uses mostly matplotlib, it is faster, less prone to errors, and you have more options available)
# for a tutorial on regular Octave-style plotting in Python, see https://matplotlib.org/stable/tutorials/introductory/pyplot.html
# for a tutorial on plotting images, see https://matplotlib.org/stable/gallery/images_contours_and_fields/image_demo.html
import matplotlib.pyplot as plt
plt.plot(c)
# e.g., for plotting images (but please check the links above for more details, this is really just a proof of concept...)
# fig, ax = plt.subplots()
# im = ax.imshow(lp.real(c))
# plt.show(c)

# b) alternatively, plot directly in Octave
lp.imagesc(abs(c));
lp.close()
# TROUBLESHOOTING PLOTS: if the plotting does not work (please note that this is highly system-dependent, so I can only give you guidelines): 
# 1. you have used somethin like 'figure(1)': depending on your system, this might cause your plot to not be visible
# 2. you have no plotting engine available:
# too see your available plotting engines for Octave, type
# lp.available_graphics_toolkits()
# to select another plotting engine, type
# lp.graphics_toolkit('TOOLKIT_NAME')

# 2.) run the wavelet-demo
lp.run('demo_wavelets')#this should give you plots, also in Octave
lp.close('all')
# 3.) audfilters-demo
#lp.help('greasy') #access the ltfat help
lp.eval("[f,fs]=greasy;")#functions with no input args do not work well when called directly
f=lp.pull('f')
L=lp.length(f)
lp.push('L', L)
# the following should be done via eval because I can not yet define function handles in Python (so g would have to be converted manually by you)
lp.eval("[g,a,fc]=audfilters(fs,L,'fractional');") 
#this needs to be executed via eval, because Python has no way of knowing that it should interpret {'realdual',g} as a cell array
fb = lp.eval("c=filterbank(f,{'realdual',g},a);")
#get the variables back to python (you could of course also leave them in Octave and calculate everything via eval)
c=lp.pull('c')
a=lp.pull('a')
fs=lp.pull('fs')
fc=lp.pull('fc')
out=lp.plotfilterbank(c,a,fc,fs,90,'audtick');
lp.close()#close the plot


# further information:
# https://oct2py.readthedocs.io/en/latest/
# https://blink1073.github.io/oct2py/source/info.html

# Python: debugging import issues
# Start your script with python -v -m my_scriptname.py and then check the output 
# to see exactly where your modules are getting imported from

