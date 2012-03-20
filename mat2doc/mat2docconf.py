
# -------------------------------------------
# Global configuration of the mat2doc system
# -------------------------------------------

import localconf

conf=ConfType()

conf.octexec=localconf.octexec
conf.matexec=localconf.matexec
conf.plotengine=localconf.plotengine
   
conf.root=localconf.userdir+'nw/ltfat/'

conf.bibfile=conf.root+'mat2doc/project'

conf.workdir=localconf.userdir+'publish/'

conf.ignorepars=['a','order','noise']

conf.copyrightplate=conf.root+'mat2doc/copyrightplate'

def mycopyrightfun(self):
    vf=file(self.root+'ltfat_version');
    v=vf.readline()
    vf.close
    
    f=file(self.copyrightplate)
    buf=f.readlines()
    f.close

    copyright=['Copyright (C) 2005-2012 Peter L. Soendergaard.\n','This file is part of LTFAT version '+v]
    copyright.extend(buf)
    
    return copyright


conf.copyright=mycopyrightfun

contentsfiles=['Contents','gabor/Contents','fourier/Contents',
               'filterbank/Contents','nonstatgab/Contents','frames/Contents',
               'sigproc/Contents','auditory/Contents',
               'demos/Contents','signals/Contents']

# ------------------------------------------
# Configuration of PHP for Sourceforge
# ------------------------------------------

php=HTMLConf()

php.basetype='html'

php.subdir='ltfathtml/'

php.indexfiles=contentsfiles

# This is the full path, used for php-including files.
php.docroot='/home/groups/l/lt/ltfat/htdocs/doc/'

# This is the usual web-server root for "<a href=" ... > tags.
php.htmlroot='/doc/'
    
php.hb='<H2>'

php.he='</H2>'

php.fext='.php'

php.head="""<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN"><html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>The Linear Time-Frequency Analysis Toolbox</title>
<link rel="stylesheet" href="/ltfat.css" type="text/css">
<script type="text/javascript"
   src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
</head>
<body>

   <?php ini_set("include_path",".:/home/groups/l/lt/ltfat/htdocs/"); ?>
   <?php require("topnav.php");?>

"""

php.foot="""<?php require("footer.php");?>
</body>
</html>"""

php.dryrun=1;


# ------------------------------------------
# Configuration of LaTeX
# ------------------------------------------

tex=TexConf()

tex.basetype='tex'

tex.subdir='toolboxref/'

tex.texfile='toolboxref'

tex.indexfiles=contentsfiles
    
tex.head="""\documentclass{amsart}
\usepackage{ae}
\usepackage{aecompl}
\usepackage[T1]{fontenc}
\usepackage[latin1]{inputenc}
\usepackage{graphicx}
\usepackage{hyperref}
\setlength\parskip{\medskipamount}
\setlength\parindent{0pt}
\makeatletter

%\usepackage{babel}
\usepackage{amssymb}
\makeatother
\\begin{document}
\\title{LTFAT Reference manual}
\\author{Peter L. S{\\o}ndergaard}

\\maketitle
\\tableofcontents{}

"""
tex.foot='\\bibliographystyle{abbrv}\n'+ \
         '\\bibliography{'+conf.bibfile+'}\n'+ \
"""\end{document}
"""

tex.dryrun=1
tex.dooutput=0

# ------------------------------------------
# Configuration of Matlab
# ------------------------------------------

mat=ConfType()

mat.basetype='mat'

mat.subdir='ltfat/'

mat.urlbase='http://ltfat.sourceforge.net/doc'
mat.urlext='php'

# ------------------------------------------
# Configuration of Verification system
# ------------------------------------------

verify=ConfType()

verify.basetype='verify'

verify.targets=['AUTHOR','TESTING','REFERENCE']

verify.notappears=['FIXME','BUG','XXL','XXX']

verify.ignore=["demo_","comp_","assert_","Contents.m","init.m"]

