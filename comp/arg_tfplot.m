function definput=arg_tfplot(definput)
  
  definput.flags.tc={'notc','tc'};
  definput.flags.plottype={'image','contour','surf','pcolor'};
  definput.flags.clim={'noclim','clim'};
  definput.flags.log={'db','dbsq','lin','linsq','linabs'};
  definput.flags.colorbar={'colorbar','nocolorbar'};
  
  definput.keyvals.fs=[];
  definput.keyvals.clim=[0,1];
  definput.keyvals.dynrange=[];
  