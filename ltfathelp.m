function op1=ltfathelp(varargin)
%LTFATHELP Help on the LTFAT toolbox
%   Usage:  ltfathelp;
%           v=ltfathelp('version');
%           mlist=ltfathelp('modules');
%
%   `ltfathelp` displays some general help on the LTFAT toolbox.
%
%   `ltfathelp('version')` returns the version number.
%
%   `ltfathelp('modules')` returns a cell array of installed modules and
%   corresponding version numbers.
%
%   See also:  ltfatstart

%   AUTHOR : Peter L. Søndergaard.  
%   TESTING: NA
%   REFERENCE: NA


  
% Verify that comp_pgauss is in path
if ~exist('comp_pgauss','file')
  disp(' ');
  disp('--- LTFAT - The Linear Time Frequency Analysis toolbox. ---');
  disp(' ')
  disp('To start the toolbox, call LTFATSTART as the first command.');
  disp(' ');
  return;
end;

bp=ltfatbasepath;

definput.keyvals.versiondata=[];
definput.keyvals.modulesdata=[];
definput.flags.mode={'general','version','modules'};

[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_general
  disp(' ');
  disp('--- LTFAT - The Linear Time Frequency Analysis toolbox. ---');
  disp(' ')

  disp(['Version ',kv.versiondata]);
  disp(' ');
  disp('Installed modules:');
  disp(' ');
  disp('Name:            Version:  Description');
  modinfo=ltfathelp('modules');
  for ii=1:length(modinfo);
    s=sprintf(' %-15s %7s  %s',modinfo{ii}.name,modinfo{ii}.version, ...
	      modinfo{ii}.description);
    disp(s);
  end;

  disp('Type "help modulename" where "modulename" is the name of one')
  disp('of the modules to see help on that module.') 

  disp(' ');
  disp('For other questions, please don''t hesitate to send an email to ltfat-help@lists.sourceforge.net.'); 
    
end;
  
if flags.do_version
  op1=kv.versiondata;
end;

if flags.do_modules
  op1={};
  for ii=1:numel(kv.modulesdata)
    
    p=kv.modulesdata{ii};
    
    % Get the first line of the help file
    [FID, MSG] = fopen ([bp,p.name,filesep,'Contents.m'],'r');
    if FID==-1
      error('Module %s does not contain a Contents.m file.',p.name);
    end;
    firstline = fgetl (FID);
    fclose(FID);
    
    
    % Load the information into the cell array.	
    op1{ii}.name=p.name;
    op1{ii}.version=p.version;
    op1{ii}.description=firstline(2:end);
  end;
end;



