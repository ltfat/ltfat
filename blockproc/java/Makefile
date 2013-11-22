# Use Make to process this file. 

JC=javac

FLAGS=-source 1.6 -target 1.6
include ../../src/ostools.mk

default:
	$(JC) $(FLAGS) net/sourceforge/ltfat/*.java net/sourceforge/ltfat/thirdparty/*.java
	jar cf blockproc.jar net/sourceforge/ltfat/*.class net/sourceforge/ltfat/thirdparty/*.class

clean:
	$(RM) net/sourceforge/ltfat/*.class
	$(RM) net/sourceforge/ltfat/thirdparty/*.class
	$(RM)  *.jar

