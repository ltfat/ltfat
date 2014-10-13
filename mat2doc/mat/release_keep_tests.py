print "Creating downloadable package"

# Remove unwanted files

s=os.path.join(conf.t.dir,'timing')
rmrf(s)
os.rmdir(s)


# Recursively remove the .git files
for root, dirs, files in os.walk(conf.t.dir, topdown=False):
    for name in files:
        if name in ['.gitattributes','.gitignore','desktop.ini']:
            os.remove(os.path.join(root, name))


# "bootstrap" the configure files
os.system("cd "+conf.t.dir+"/src; ./bootstrap")

s=os.path.join(conf.t.dir,'src','autom4te.cache')
rmrf(s)
os.rmdir(s)

# Compile the Java classes
os.system("cd "+conf.t.dir+"/blockproc/java; make")
os.system("cd "+conf.t.dir+"/blockproc/java; make classclean")
