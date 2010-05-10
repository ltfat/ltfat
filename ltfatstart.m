function ltfatstart()
%LTFATSTART   Start the LTFAT toolbox.
%   Usage:  ltfatstart;
%
%   LTFATSTART starts the LTFAT toolbox. This command must be run
%   before using any of the functions in the toolbox. If you issue a
%   CLEAR ALL then you must call LTFATSTART again.
%
%   You may manually edit this function to configure the toolbox to your
%   taste. This includes default settings for the plotting functions.
%
%   See also:  ltfathelp

%   AUTHOR : Peter Soendergaard.  
%   TESTING: NA
%   REFERENCE: NA

% ------- Settings: ----------------------------------

% --- general settings ---
% Print the banner at startup?
printbanner=1;

% --- Settings for iterative algorithms ---

% Problem sizes below this value will be handled by a direct approach.
itercrossover=200;

% Maximum number of iterations to do:
itermaxit=200;

% The precision to stop at:
itertol=1e-10;

% --- Settings for plots ---

% Approximate resolution along time-axis.
xres=800;

% Ratio of frequency resolution versus time resolution
displayratio=3/4;


% ----------------------------------------------------
% -------   do not edit below this line   ------------
% ----------------------------------------------------

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('ltfatstart');
% Kill the function name from the path.
basepath=basepath(1:end-13);

% add the base path
addpath(basepath);

bp=[basepath,filesep];

% Load the version number
[FID, MSG] = fopen ([bp,'ltfat_version'],'r');
if FID == -1
    error(MSG);
else
    ltfat_version = fgetl (FID);
    fclose(FID);
end

% Create and load the information.
global TF_CONF;
TF_CONF.basepath=bp;
TF_CONF.ltfat_version=ltfat_version;

% -----------  install the modules -----------------

modules={};
nplug=0;

% List all files in base directory
d=dir(basepath);

for ii=1:length(d)
  
  % We only look for directories
  if ~d(ii).isdir
    continue;
  end;
  
  % Skip the default directories . and ..
  if (d(ii).name(1)=='.')
    continue;
  end;
  
  % Skip directories without an init file
  name=d(ii).name;
  % CThe file is a directory and it does not start with '.' This could
  % be a module
  if ~exist([bp,name,filesep,name,'init.m'],'file')
    continue
  end;
    
  % Now we know that we have found a module
  
  % Set 'status' to zero if the module forgets to define it.
  status=0;
  
  module_version=ltfat_version;
     
  % Add the module dir to the path
  addpath([bp,name])  
  
  % Execute the init file to see if the status is set.
  eval([name,'init']);
  if status>0
    if status==1
      nplug=nplug+1;
      modules{nplug}.name=name;
      modules{nplug}.version=module_version;
    end;
  else
    % Something failed, restore the path
    rmpath([bp,name]);
  end;
end;

TF_CONF.modules=modules;

% Check if Octave was called using 'silent'
%if isoctave
%  args=argv;
%  for ii=1:numel(args)
%    s=lower(args{ii});
%    if strcmp(s,'--silent') || strcmp(s,'-q')
%      printbanner=0;
%    end;
%  end;
%end;

if printbanner
  s=which('comp_pgauss');
  if isempty(s)
    error('comp_pgauss not found, something is wrong.')
  end;
  
  if strcmp(s(end-1:end),'.m')
    backend = 'LTFAT is using the script language backend.';
  else
    if isoctave
      backend = 'LTFAT is using the C++ Octave backend.';
    else
      backend = 'LTFAT is using the MEX backend.';
    end;
  end;
  
  banner = sprintf(['LTFAT version %s. Copyright 2010 Peter L. Soendergaard. ' ...
                    'For help, please type "ltfathelp". %s'], ...
                   ltfat_version,backend);
  
  disp(banner);
end;
