#!/usr/bin/python

import sys,os
cwd=os.getcwd()+os.sep

# ------- Configuration parameters -------------

# Configure HTML placement at remote server
host='soender,ltfat@web.sourceforge.net'
www='/home/project-web/ltfat/htdocs//'

tbwww='/home/peter/nw/ltfatwww'

# ------- do not edit below this line ----------

# Import the localconf file, it should be placed in the same directory
# as this function is called from.

sys.path.append(cwd)
import localconf

sys.path.append(localconf.mat2docdir)
import mat2doc
import printdoc

# Safely create directories if they are missing
printdoc.safe_mkdir(localconf.outputdir)
printdoc.safe_mkdir(localconf.tmpdir)

# This is run if the file is called from the command line
project=sys.argv[1]
todo=sys.argv[2]
redomode='auto'
if len(sys.argv)>3:
    redomode=sys.argv[3]


sys.path.append(localconf.projects[project]+'mat2doc')
from mat2docconf import *

filesdir=localconf.outputdir+project+'-files'+os.sep
printdoc.safe_mkdir(filesdir)

def createcompressedfile(srcdir,fname,target):
    # The following targets are supported:
    # unix - Unix lineendings, UTF-8 encoding and gzip compression
    # win  - Windows lineendings, 8859 enconding and zip compression
    # mac  - Unix lineendings, 8859 encoding and zip compression
    tmpworkdir=localconf.tmpdir+project
    fname=filesdir+fname
    srcdir=localconf.outputdir+srcdir+os.sep

    printdoc.safe_mkdir(tmpworkdir)
    printdoc.rmrf(tmpworkdir)

    # Create zip and tgz archives.  In order to get the right path in
    # the compressed files, copy the source to the temporary directory
    # with the project
    printdoc.rmrf(tmpworkdir)
    os.system('cp -R '+srcdir+'* '+tmpworkdir)
    
    if target=='unix':
        os.system('cd '+localconf.tmpdir+'; tar zcvf '+fname+'.tgz '+project)
    
    if target=='win':
        os.system('rm '+fname+'.zip')
        printdoc.unix2dos(tmpworkdir)
        printdoc.convertencoding(tmpworkdir,'ISO-8859-1')
        os.system('cd '+localconf.tmpdir+'; zip -r '+fname+'.zip '+project)

    if target=='mac':
        os.system('rm '+fname+'.zip')
        printdoc.convertencoding(tmpworkdir,'ISO-8859-1')
        os.system('cd '+localconf.tmpdir+'; zip -r '+fname+'.zip '+project)
    
def createbinaryfile(filename,ext,target):
    tmpworkdir=localconf.outputdir+project+'-mat'        
    bdir=localconf.outputdir+project+'-'+ext

    printdoc.safe_mkdir(bdir)
    printdoc.rmrf(bdir)
    
    s='cp -r '+tmpworkdir+'/* '+bdir
    os.system(s)

    s='cp -r '+localconf.outputdir+project+'-'+ext+'-addon/* '+bdir
    os.system(s)

    createcompressedfile(project+'-'+ext,project+'-'+versionstring+'-'+ext,target)

def runcommand(todo,redomode='auto'):
    # When editing these commands, some variables are already defined
    #
    #   versionstring - This is the string from mat2docconf
    #   project       - This is the name of the project from calling publish.py

    print 'PUBLISH '+todo+' '+redomode

    # Simple commands to make the others easier to write
    if todo=='svnmat':
        pass

    if todo=='gitstagemat':
        printdoc.git_stageexport_mat(project)
        printdoc.printdoc(project,'mat')

    if todo=='gitrepomat':
        printdoc.git_repoexport_mat(project)
        printdoc.printdoc(project,'mat')


    # Release for other developers to download
    if 'develmat' in todo:
        runcommand('gitstagemat')

        createcompressedfile(project+'-mat',project+'-devel-'+versionstring,'unix')
        createcompressedfile(project+'-mat',project+'-devel-'+versionstring,'win')

    # Release for users to download
    if 'releasemat' in todo:
        runcommand('gitrepomat')

        matdir=localconf.outputdir+project+'-mat'+os.sep

        # Remove unwanted files
        os.system('rm -rf '+matdir+'testing')
        os.system('rm -rf '+matdir+'reference')
        os.system('rm -rf '+matdir+'timing')

        createcompressedfile(project+'-mat',project+'-'+versionstring,'unix')
        createcompressedfile(project+'-mat',project+'-'+versionstring,'win')

    if 'tex'==todo:
        printdoc.printdoc(project,'tex')

    if 'texmake'==todo:
        os.system('cd '+project['tex']+'; make')

        printdoc.printdoc(project,'tex')

    if 'texupload'==todo:
        texdir=localconf.outputdir+project+'-tex'+os.sep
        s='rsync -av '+texdir+'ltfat.pdf '+host+':'+www+'doc/'
        os.system(s)

    if todo=='php':
        phpdir=localconf.outputdir+project+'-php'+os.sep
        printdoc.printdoc(project,'php')
        s='rsync -av '+phpdir+' '+host+':'+www+'doc/'
        os.system(s)    

    if todo=='phplocal' in todo:
        printdoc.printdoc(project,'phplocal')

    if todo=='fullrelease':
        runcommand('releasemat',redomode)
        runcommand('binary')
        runcommand('php',redomode)
        runcommand('tex',redomode)
        runcommand('texupload')
    
    if todo=='wavephp':
        wavedir=localconf.projects['ltfatwave']
        printdoc.assert_git_on_branch(wavedir,'wavelets')
        printdoc.printdoc('ltfatwave','php')

    if todo=='wavestagemat': 
        printdoc.git_stageexport_mat('ltfatwave')
        printdoc.printdoc('ltfatwave','mat')

    if 'stagewww'==todo:
        publishwww=localconf.outputdir+'ltfatwww/'
        printdoc.git_stageexport(tbwww,publishwww)
        os.system('rsync -av '+publishwww+' '+host+':'+www);

    if 'releasewww'==todo:
        publishwww=localconf.outputdir+'ltfatwww/'
        printdoc.git_repoexport(tbwww,'master','ltfatwww',localconf.outputdir)
        os.system('rsync -av '+publishwww+' '+host+':'+www);

    if 'binary'==todo:
        # You must run this command right after the "releasemat" or
        # "develmat" commands as this will create a correct "-mat"
        # directory

        createbinaryfile(project+'-'+versionstring,'win64','win')
        createbinaryfile(project+'-'+versionstring,'mac','mac')


    #if 'upload' in todo:
    #    ddir=localconf.outputdir+'ltfat_sourceforge/ltfat/'
    #    os.system('rsync -av '+ddir+
    #              ' soender,ltfat@frs.sourceforge.net:/home/frs/project/l/lt/ltfat/ltfat/')

    if 'notesmake'==todo:
        notes=notes.getnotenumbers(notesdir)

        notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

        for notenumber in notes:
            print 'Trying to make LTFAT note '+notenumber
            os.system('cd '+notesdir+notenumber+'; make')

    if 'notestexclean'==todo:
        notes=notes.getnotenumbers(notesdir)

        notes = filter(lambda x: (os.path.exists(notesdir+x+'/Makefile')), notes)

        for notenumber in notes:
            os.system('cd '+notesdir+notenumber+'; make texclean')

    if 'noteshtml'==todo:

        printdoc.printnoteshtml('ltfatnote',notesdir,notehtml)

        os.system('rsync -av '+notehtml+' '+host+':'+noteswww);




# Excute runcommand, this is where the main stuff happens    
runcommand(todo,redomode)

