print "Creating downloadable package"

# Remove unwanted files
s=os.path.join(conf.t.dir,'testing')
rmrf(s)
os.rmdir(s)

s=os.path.join(conf.t.dir,'timing')
rmrf(s)
os.rmdir(s)

s=os.path.join(conf.t.dir,'reference')
rmrf(s)
os.rmdir(s)

os.remove(os.path.join(conf.t.dir,'.gitattributes'))
os.remove(os.path.join(conf.t.dir,'.gitignore'))



