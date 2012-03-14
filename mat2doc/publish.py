#!/usr/bin/python

import sys,os

cwd=os.getcwd()+'/'

sys.path.append(cwd)

from localconf import *
from datesuffix import *


tbpath  =basepath+'ltfat/'

# Configure HTML placement at remote server
host='soender,ltfat@web.sourceforge.net'
www='/home/groups/l/lt/ltfat/htdocs/'

# do not edit below this line

tbwww   =basepath+'ltfatwww/'
m2dfile =tbpath+'/mat2doc/mat2docconf.py'

publishmat=cwd+'ltfat/'
notesdir=basepath+'notes/'
notehtml=cwd+'noteshtml/'
noteswww=www+'notes/'

f=file(tbpath+'ltfat_version')
versionstring=f.read()[:-1]
f.close()

sys.path.append(basepath+'mat2doc/')

import printdoc, notes

todo=sys.argv[1]

if 'verify' in todo or todo==[]:
    printdoc.printdoc([m2dfile,'verify'])

if 'stagemat' in todo:
    printdoc.git_stageexport(tbpath,publishmat)
    printdoc.print_mat(m2dfile,publishmat)

    fname=cwd+'ltfat-'+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'ltfat')
    os.system('tar zcvf '+fname+'.tgz ltfat/')

    printdoc.unix2dos(cwd+'ltfat')
    os.system('zip -r '+fname+'.zip ltfat/')


if 'releasemat' in todo:
    printdoc.git_repoexport(tbpath,'master','ltfat',cwd)
    printdoc.print_mat(m2dfile,publishmat)
    
    # Remove unwanted files
    os.system('rm -rf '+publishmat+'testing')
    os.system('rm -rf '+publishmat+'reference')
    os.system('rm -rf '+publishmat+'timing')

    fname=cwd+'ltfat-'+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'ltfat')
    os.system('tar zcvf '+fname+'.tgz ltfat/')

    printdoc.unix2dos(cwd+'ltfat')
    os.system('zip -r '+fname+'.zip ltfat/')

if 'develmat' in todo:
    printdoc.git_repoexport(tbpath,'master','ltfat',cwd)

    fname=cwd+'ltfat-devel-'+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'ltfat')
    os.system('tar zcvf '+fname+'.tgz ltfat/')

    printdoc.unix2dos(cwd+'ltfat')
    os.system('zip -r '+fname+'.zip ltfat/')

if 'releasebranch' in todo:
    bname=sys.argv[2]
    printdoc.git_repoexport(tbpath,bname,'ltfat',cwd)
    printdoc.print_mat(m2dfile,publishmat)
    
    # Remove unwanted files
    os.system('rm -rf '+publishmat+'testing')
    os.system('rm -rf '+publishmat+'reference')
    os.system('rm -rf '+publishmat+'timing')

    fname=cwd+'ltfat-'+bname+'-'+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'ltfat')
    os.system('tar zcvf '+fname+'.tgz ltfat/')

    printdoc.unix2dos(cwd+'ltfat')
    os.system('zip -r '+fname+'.zip ltfat/')

if 'releasetestbranch' in todo:
    bname=sys.argv[2]
    printdoc.git_repoexport(tbpath,bname,'ltfat',cwd)
    printdoc.print_mat(m2dfile,publishmat)
    
    fname=cwd+'ltfat-'+bname+'-'+versionstring
    os.system('rm '+fname+'.zip')

    printdoc.dos2unix(cwd+'ltfat')
    os.system('tar zcvf '+fname+'.tgz ltfat/')

    printdoc.unix2dos(cwd+'ltfat')
    os.system('zip -r '+fname+'.zip ltfat/')

    
if 'pdf' in todo:
    printdoc.printdoc([m2dfile,'tex'])
    os.system('cd toolboxref; pdflatex toolboxref.tex')
    os.system('cd toolboxref; bibtex toolboxref')
    os.system('cd toolboxref; pdflatex toolboxref.tex')
    
    os.system('scp toolboxref/toolboxref.pdf '+host+':'+www+'doc/ltfat.pdf')

if 'stagewww' in todo:
    publishwww=cwd+'ltfatwww/'
    printdoc.autostage(tbwww)
    printdoc.git_stageexport(tbwww,publishwww)
    os.system('cp ltfat-devel-'+versionstring+'.zip '+publishwww+'/prerelease/')    
    os.system('cp ltfat-devel-'+versionstring+'.tgz '+publishwww+'/prerelease/')
    os.system('cp ltfat-devel-'+versionstring+'-win32.zip '+publishwww+'/prerelease/')

    os.system('rsync -av '+publishwww+' '+host+':'+www);

if 'releasewww' in todo:
    publishwww=cwd+'ltfatwww/'
    printdoc.git_repoexport(tbwww,'master','ltfatwww',cwd)
    os.system('rsync -av '+publishwww+' '+host+':'+www);

if 'binary' in todo:
    # Build windows binary
    fname=cwd+'ltfat-'+versionstring+'-win32'
    os.system('rm '+fname+'.zip')

    bdir=cwd+'buildbinary/ltfat'
    printdoc.rmrf(bdir)

    printdoc.unix2dos(cwd+'ltfat')
    
    os.system('cp -r '+cwd+'ltfat/* '+bdir)
    os.system('cp -r '+cwd+'ltfat-win32-addon/* '+bdir)

    s='cd '+cwd+'buildbinary; zip -r '+fname+'.zip ltfat/'
    print s
    os.system(s)

    # Build Mac binary
    fname=cwd+'ltfat-'+versionstring+'-mac'
    os.system('rm '+fname+'.zip')

    bdir=cwd+'buildbinary/ltfat'
    printdoc.rmrf(bdir)

    printdoc.unix2dos(cwd+'ltfat')
    
    os.system('cp -r '+cwd+'ltfat/* '+bdir)
    os.system('cp -r '+cwd+'ltfat-mac-addon/* '+bdir)

    s='cd '+cwd+'buildbinary; zip -r '+fname+'.zip ltfat/'
    print s
    os.system(s)



    

if 'php' in todo:
    printdoc.printdoc([m2dfile,'php'])
    os.system('cp '+tbwww+'doc/index.php ltfathtml/')
    s='rsync -av ltfathtml/ '+host+':'+www+'doc/'
    os.system(s)

if 'phpupload' in todo:
    s='rsync -av ltfathtml/ '+host+':'+www+'doc/'
    os.system(s)    

if 'upload' in todo:
    ddir=cwd+'ltfat_sourceforge/ltfat/'
    os.system('rsync -av '+ddir+
              ' soender,ltfat@frs.sourceforge.net:/home/frs/project/l/lt/ltfat/ltfat/')

if 'notesmake' in todo:
    notes=notes.getnotenumbers(notesdir)

    notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

    for notenumber in notes:
        print 'Trying to make LTFAT note '+notenumber
        os.system('cd '+notesdir+notenumber+'; make')

if 'notestexclean' in todo:
    notes=notes.getnotenumbers(notesdir)

    notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

    for notenumber in notes:
        os.system('cd '+notesdir+notenumber+'; make texclean')

if 'noteshtml' in todo:

    printdoc.printnoteshtml('ltfatnote',notesdir,notehtml)
        
    os.system('rsync -av '+notehtml+' '+host+':'+noteswww);




