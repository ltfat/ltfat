# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

# When writing this file, certain variables are already defined:
#
#   self.root points to the project directory


import localconf

conf=ConfType()

conf.urlbase='http://ltfat.sourceforge.net/doc/'

def mycopyrightfun(self):
    vf=file(self.root+'ltfat_version');
    v=vf.readline()
    vf.close
    
    f=file(self.root+'mat2doc/copyrightplate')
    buf=f.readlines()
    f.close

    copyright=['Copyright (C) 2005-2012 Peter L. Soendergaard.\n','This file is part of LTFAT version '+v]
    copyright.extend(buf)
    
    return copyright

conf.copyright=mycopyrightfun

contentsfiles=['Contents','gabor/Contents','fourier/Contents',
               'filterbank/Contents','nonstatgab/Contents',
               'frames/Contents',
               'sigproc/Contents','auditory/Contents',
               'demos/Contents','signals/Contents']

# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=PhpConf()

php.indexfiles=contentsfiles

# This is the full path, used for php-including files.
php.docroot='/home/project-web/ltfat/htdocs/doc/'

# This is the usual web-server root for "<a href=" ... > tags.
php.fext='.php'

php.head="""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN"><html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>The Linear Time-Frequency Analysis Toolbox</title>
<meta> 
<link rel="stylesheet" href="/ltfat.css" type="text/css">
<script type="text/javascript"
   src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
</head>
<body>

   <?php ini_set("include_path",".:/home/project-web/ltfat/htdocs/"); ?>
   <?php include("topnav.php");?>

"""

php.foot="""<?php include("footer.php");?>
</body>
</html>"""


# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.indexfiles=contentsfiles
    
# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=MatConf()

# ------------------------------------------
# Configuration of Verification system
# ------------------------------------------

verify=ConfType()

verify.basetype='verify'

verify.targets=['AUTHOR','TESTING','REFERENCE']

verify.notappears=['FIXME','BUG','XXL','XXX']

verify.ignore=["demo_","comp_","assert_","Contents.m","init.m"]



