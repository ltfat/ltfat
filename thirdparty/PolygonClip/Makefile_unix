# This line gets the MEX extension for the platform. If the shell
# command 'mexext' is not available on your Matlab installation,
# please execute the 'mexext' command within Matlab, and put the
# answer into the line below (remember to uncomment).
EXT = $(shell mexext)
#EXT = $(glnx86)

mex:
	mex -I../GPC PolygonClip.c ../GPC/gpc.c

clean:
	rm PolygonClip.$(EXT)