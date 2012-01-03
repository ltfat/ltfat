function op1=ltfathelp(varargin)
%LTFATHELP Help on the LTFAT toolbox.
%   Usage:  ltfathelp;
%           v=ltfathelp('version');
%           mlist=ltfathelp('modules');
%
%   LTFATHELP displays some general help on the LTFAT toolbox.
%
%   LTFATHELP('version') returns the version number.
%
%   LTFATHELP('modules') returns a cell array of installed modules and
%   corresponding version numbers.
%
%   LTFATHELP('authors') lists the name of the authors.
%
%   See also:  ltfatstart

%   AUTHOR : Peter Soendergaard.  
%   TESTING: NA
%   REFERENCE: NA


  
% Verify that ltfatarghelper is in path
if ~exist('ltfatarghelper','file')
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
definput.flags.mode={'general','version','modules','authors'};

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

  disp(' ')
  if isoctave
    disp('Type ltfathelp("modulename") where "modulename" is the name of one');
    disp('of the modules to see help on that module.');
    
  else
    disp('Type "help modulename" where "modulename" is the name of one')
    disp('of the modules to see help on that module.') 

  end; 
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

if flags.do_authors
  disp('Peter L. Soendergaard, Centre for Applied Hearing Research,');
  disp('                       Technical University of Denmark.');
  disp(' ');
  disp('Bruno Torresani, LABORATOIRE D''ANALYSE, TOPOLOGIE ET PROBABILITES');
  disp('                 Universite de Provence, Marseille, France');
  disp(' ');
  disp('Peter Balazs, Acoustics Research Institute');
  disp('              Austrian Academy of Sciences, Vienna, Austria');
  disp(' ');
  disp('Florent Jaillet, IRISA');
  disp('              Rennes, France');
end;


if 0

  bp=TF_CONF.basepath;

  switch(lower(q))
    case {'v','version'}
      if nargin==1
	op1=TF_CONF.ltfat_version;
      else
	for ii=1:length(TF_CONF.modules);
      % strcmpi does not exist in older versions of Octave,
      % therefore use strcmp(lower(...
	  if strcmp(lower(modulename),TF_CONF.modules{ii}.name)
	    found=1;
	    op1=TF_CONF.modules{ii}.version;
	    return
	  end;
	end;

	% If we get here, it means that no module matched 'modulename'.
	error('Unknown module name.');
				   
      end;
      
    otherwise
      found=0;
      for ii=1:length(TF_CONF.modules);
	if strcmp(lower(q),TF_CONF.modules{ii}.name)
	  found=1;
	  p=TF_CONF.modules{ii};

	  if isoctave
	    s=pwd;
	    cd([bp,p.name]);
	    help Contents
	    cd(s);
	  else
	    help(lower(q))
	  end;
	end;
      end;  
      if ~found
	error('Unknown command or module name.');
      end;
  end;
end;

%OLDFORMAT
