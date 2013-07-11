print "Uploading the php files"

s='rsync -av '+conf.t.dir+'/ '+conf.g.username+',ltfat@web.sourceforge.net:/home/project-web/ltfat/htdocs/doc/'
print '   '+s
os.system(s)    
