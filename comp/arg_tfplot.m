function definput=arg_tfplot(definput)
  
  definput.flags.tc={'notc','tc'};
  definput.flags.plottype={'image','contour','surf','pcolor'};  
  definput.flags.log={'db','dbsq','lin','linsq','linabs'};
  definput.flags.colorbar={'colorbar','nocolorbar'};
  definput.flags.display={'display','nodisplay'};

  definput.keyvals.fontsize=14;  
  definput.keyvals.fs=[];
  definput.keyvals.clim=[];
  definput.keyvals.dynrange=[];  
  definput.keyvals.colormap=ltfat_inferno(); 
