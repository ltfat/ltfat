#!/usr/bin/python

import sys,os
cwd=os.getcwd()+os.sep

# ------- Configuration parameters -------------

projectname='ltfat'

# Configure HTML placement at remote server
host='soender,ltfat@web.sourceforge.net'
www='/home/project-web/ltfat/htdocs//'

tbwww='/home/peter/nw/ltfatwww'

# ------- do not edit below this line ----------

# Import the localconf file, it should be place in the same directory
# as this function is called from.

sys.path.append(cwd)

import localconf

sys.path.append(localconf.mat2docdir)

# Get the data from localconf
projectdir=localconf.projects[projectname]
conffile=projectdir+'mat2doc/mat2docconf.py'
filesdir=localconf.outputdir+projectname+'-files'+os.sep

f=file(projectdir+projectname+'_version')
versionstring=f.read()[:-1]
f.close()

import printdoc

def runcommand(todo,redomode='auto'):

    print 'PUBLISH '+todo+' '+redomode

    # Release for other developers to download
    if 'develmat' in todo:
        printdoc.git_stageexport_mat(projectname)
        printdoc.printdoc(projectname,'mat')

        fname=filesdir+projectname+'-devel-'+versionstring
        os.system('rm '+fname+'.zip')

        # Create the Unix src package
        os.system('tar zcvf '+fname+'.tgz '+projectname+os.sep)

        # Create the Windows src package
        os.system('rm '+fname+'.zip')
        printdoc.unix2dos(filesdir+projectname)
        os.system('zip -r '+fname+'.zip '+projectname+os.sep)

    # Release for users to download
    if 'releasemat' in todo:
        printdoc.git_repoexport_mat(projectname)
        matdir=localconf.outputdir+projectname+'-mat'+os.sep

        # Remove unwanted files
        os.system('rm -rf '+matdir+'testing')
        os.system('rm -rf '+matdir+'reference')
        os.system('rm -rf '+matdir+'timing')

        printdoc.printdoc(projectname,'mat')

        fname=filesdir+projectname+'-'+versionstring
        printdoc.safe_mkdir(filesdir)

        # Create the binary packages In order to get the right path in
        # the compressed files, copy the source to the temporary
        # directory with the projectname
        os.system('cp -R '+matdir+' '+localconf.tmpdir+projectname)

        # Create the Unix src package
        os.system('cd '+localconf.tmpdir+'; tar zcvf '+fname+'.tgz '+projectname)

        # Create the Windows src package
        os.system('rm '+fname+'.zip')
        printdoc.unix2dos(filesdir+projectname)
        os.system('cd '+localconf.tmpdir+'; zip -r '+fname+'.zip '+projectname)

    if 'tex'==todo:
        printdoc.printdoc(projectname,'tex')

    if 'texmake'==todo:
        os.system('cd '+project['tex']+'; make')

        printdoc.printdoc(projectname,'tex')

    if 'texupload'==todo:
        texdir=localconf.outputdir+projectname+'-tex'+os.sep
        s='rsync -av '+texdir+'ltfat.pdf '+host+':'+www+'doc/'
        os.system(s)

    if todo=='php':
        phpdir=localconf.outputdir+projectname+'-php'+os.sep
        printdoc.printdoc(projectname,'php')
        s='rsync -av '+phpdir+' '+host+':'+www+'doc/'
        os.system(s)    

    if todo=='phplocal' in todo:
        printdoc.printdoc(projectname,'phplocal')

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
        publishwww=cwd+'ltfatwww/'
        printdoc.git_stageexport(tbwww,publishwww)
        #os.system('cp ltfat-devel-'+versionstring+'.zip '+publishwww+'/prerelease/')    
        #os.system('cp ltfat-devel-'+versionstring+'.tgz '+publishwww+'/prerelease/')
        #os.system('cp ltfat-devel-'+versionstring+'-win32.zip '+publishwww+'/prerelease/')

        os.system('rsync -av '+publishwww+' '+host+':'+www);

    if 'releasewww'==todo:
        publishwww=cwd+'ltfatwww/'
        printdoc.git_repoexport(tbwww,'master','ltfatwww',cwd)
        os.system('rsync -av '+publishwww+' '+host+':'+www);

    if 'binary'==todo:
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


    #if 'upload' in todo:
    #    ddir=cwd+'ltfat_sourceforge/ltfat/'
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


# This is run if the file is called from the command line
todo=sys.argv[1]
redomode='auto'
if len(sys.argv)>2:
    redomode=sys.argv[2]
    
runcommand(todo,redomode)


