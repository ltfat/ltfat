print("Creating downloadable package")

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

# Recursively remove the .git files
for root, dirs, files in os.walk(conf.t.dir, topdown=False):
    for name in files:
        if name in ['.gitattributes','.gitignore','desktop.ini']:
            os.remove(os.path.join(root, name))


# Compile the Java classes
os.system("cd "+conf.t.dir+"/blockproc/java; make")
os.system("cd "+conf.t.dir+"/blockproc/java; make classclean")
