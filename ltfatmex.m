function ltfatmex(optclean)
%LTFATMEX   Compile Mex/Oct interfaces
%   Usage:  ltfatmex;
%           ltfatmex('clean');
%
%   LTFATMEX compiles certain functions written in C in order to speed
%   up the execution of the toolbox.
%
%   LTFATMEX('clean') removes the compiled functions.

%   AUTHOR : Peter Soendergaard.
%   TESTING: NA
%   REFERENCE: NA

error(nargchk(0,1,nargin));

global TF_CONF;

% Verify that TF_CONF has been initialized
if numel(TF_CONF)==0
    disp('');
    disp('--- LTFAT - The Linear Time Frequency Analysis toolbox. ---');
    disp('')
    disp('To start the toolbox, call LTFATSTART as the first command.');
    disp('');
    return;
end;

bp=TF_CONF.basepath;

% Remember the current directory.
curdir=pwd;

if (nargin>0) && strcmpi(optclean,'clean')
    
    if isoctave
        cd([bp,'oct']);
        
        deletefiles('*.oct');
        deletefiles('*.o');
        
        cd(curdir);
        
    else
        
        % Delete files permanently (bypass trashcan on Windows/Mac
        % but remember the old state
        oldstate = recycle('off');
        
        cd([bp,'mex']);
        
        deletefiles(['*.',mexext]);
        
        recycle(oldstate);

        cd(curdir);

    end;
    
    return;
    
end;

if ispc
  % FIXME Need to differentiate between Matlab and Octave.
  s=[bp,'lib',filesep,'libltfat-nomem.lib'];
else
  s=[bp,'lib',filesep,'libltfat-nomem.a'];
end;  

if ~exist(s,'file')
  error(['The LTFAT C library cannot be found. Please compile as described ' ...
         'in the INSTALL file, or download a binary release.']);
end;

if ispc && ~isoctave
  
  if ~exist([bp,'mex\libfftw3-3.dll'],'file')
    error(['Please download the FFTW binary release for Windows as ' ...
           'described in the INSTALL file.']);
  end;
  
  % Create the .lib file for the Matlab FFTW3 lib.
  if ~exist([bp,'mex\libfftw3-3.lib'],'file')
    system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
            bp,'mex\LIBFFTW3-3.DLL"']);
  end;
  
  if ~exist([bp,'mex\libfftw3f-3.dll'],'file')
    error(['Please download the FFTW binary release for Windows as ' ...
           'described in the INSTALL file.']);
  end;
  
  if ~exist([bp,'mex\libfftw3f-3.lib'],'file')
    system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
            bp,'mex\libfftw3f-3.dll"']);
  end;
  
end;

if isoctave
  cd([bp,'oct']);
  
  ext='oct';
  
  % Get the list of files.
  L=dir('*.cc');
  
  endchar=2;
else
  
  cd([bp,'mex']);
  
  ext=mexext;
  
  % Get the list of files.
  L=dir('comp_*.c');
  
  endchar=1;
end;

for ii=1:numel(L)
  filename = L(ii).name;
  objname  = [filename(1:end-endchar),ext];
  objdirinfo = dir(objname);
  
  % Make-like behaviour: build only the files where the src file is
  % newer than the object file, or the object file is missing.
  
  if isempty(objdirinfo) || (objdirinfo.datenum<L(ii).datenum)
    
    fprintf('Compiling %s\n',filename);
    
    if isoctave
      mkoctfile('-c','oct-memalloc.c');
      mkoctfile('-I../thirdparty',...
                '-I.','-I../src','-L../src','-L../lib',...
                filename,'oct-memalloc.o',...
                '-lltfat-nomem',...
                '-lfftw3',...
                '-lfftw3f',...
                '-llapack',...
                '-lblas');
      
    else
      if ispc
        mex('-c','mex-memalloc.c');
        mex('-I../thirdparty',...
            '-I.','-I../src','-L../src',...
            filename,'mex-memalloc.obj',...
            '../lib/libltfat-nomem.lib',...
            '../mex/LIBFFTW3-3.LIB',...
            '../mex/LIBFFTW3F-3.LIB',...
            [matlabroot '\extern\lib\win32\lcc\libmwlapack.lib'],...
            [matlabroot '\extern\lib\win32\lcc\libmwblas.lib']);
        
        
      else
        mex('-c','mex-memalloc.c');
        mex('-I../thirdparty',...
            '-I.','-I../src','-L../src','-L../lib',...
            filename,'mex-memalloc.o',...
            '-lltfat-nomem',...
            '-lfftw3',...
            '-lfftw3f',...
            '-llapack',...
            '-lblas');
        
      end;
      
    end;
    
    
  end;        
  
end;


if ~isoctave
  if exist([bp,filesep,'thirdparty',filesep,'GPC',filesep,'gpc.c'], ...
           'file')
    % Compile the PolygonClip interface to GPC for use with mulaclab
    cd([bp,filesep,'thirdparty',filesep,'PolygonClip']);
    mex('-I../GPC','PolygonClip.c','../GPC/gpc.c');
    
  end;
end;


% Jump back to the original directory.
cd(curdir);


function deletefiles(files)

L=dir(files);

for ii=1:numel(L)
    delete(L(ii).name);
end;

