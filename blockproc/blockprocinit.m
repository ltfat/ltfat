status=1;

% Add entry to the dynamic classpath if JVM is present.
if ~isempty(which('javaaddpath')) 
   try
      javaaddpath([basepath,filesep,'blockproc',filesep,'java',filesep]);
   catch 
       % Use lasterr for Octave compatibility
       err=lasterr;
       warning('%s: JVM support not present.',upper(mfilename));
   end
else
   warning('%s: Java toolbox not present.',upper(mfilename));
end