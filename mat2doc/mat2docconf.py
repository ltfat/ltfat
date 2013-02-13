# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

# When writing this file, the following variables are available
#
#   projectname  - The name of your project, OK, you already know this
#   projectdir   - The directory where your project is located
#   outputdir    - The output directory

from mat2doc import *

# Define versionstring
f=file(projectdir+'ltfat_version')
versionstring=f.read()[:-1]
f.close

# Define copyright    
f=file(projectdir+'mat2doc/copyrightplate')
buf=f.readlines()
f.close

copyright=[u'Copyright (C) 2005-2013 Peter L. S\xf8ndergaard.\n',
           u'This file is part of LTFAT version '+versionstring+'\n']
copyright.extend(buf)
    

conf=ConfType()

conf.copyright=copyright

contentsfiles=['Contents','gabor/Contents','fourier/Contents',
               'filterbank/Contents','nonstatgab/Contents',
               'frames/Contents',
               'sigproc/Contents','auditory/Contents','wavelets/Contents',
               'demos/Contents','signals/Contents']

# The urlbase in the targets must always be an absolute path, and it
# must end in a slash

# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=PhpConf()

php.indexfiles=contentsfiles
php.includedir='../include/'
php.urlbase='/doc/'
php.codedir=outputdir+'ltfat-mat'+os.sep

# ------------------------------------------
# Local php
# ------------------------------------------

phplocal=PhpConf()
phplocal.indexfiles=contentsfiles
phplocal.includedir='../include/'
phplocal.urlbase='/doc/'
phplocal.codedir=outputdir+'ltfat-mat'+os.sep

# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

# No demos
texcontentsfiles=['Contents','gabor/Contents','fourier/Contents',
               'filterbank/Contents','nonstatgab/Contents',
               'frames/Contents',
               'sigproc/Contents','auditory/Contents',
               'signals/Contents']


tex.indexfiles=contentsfiles
tex.urlbase='http://ltfat.sourceforge.net/doc/'
tex.codedir=outputdir+'ltfat-mat'+os.sep
    
# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=MatConf()
mat.urlbase='http://ltfat.sourceforge.net/doc/'

# ------------------------------------------
# Configuration of Verification system
# ------------------------------------------

verify=ConfType()
verify.basetype='verify'
verify.sourcedir=outputdir+'ltfat-mat'+os.sep

verify.targets=['AUTHOR','TESTING']
verify.notappears=['FIXME','BUG','XXL','XXX']

verify.ignore=["comp_","assert_","Contents.m","init.m"]



