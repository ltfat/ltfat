function ltfatmex(varargin)
%LTFATMEX   Compile Mex/Oct interfaces
%   Usage:  ltfatmex;
%           ltfatmex(...);
%
%   LTFATMEX compiles the C backend in order to speed up the execution of
%   the toolbox. The C backend is linked to Matlab and Octave through mex
%   and Octave C++ interfaces.
%
%   The action of LTFATMEX is determined by one of the following flags:
%
%-     'compile' - Compile stuff. This is the default.
%
%-     'clean'   - Removes the compiled functions.
%
%-     'test'    - Run some small tests that verify that the compiled
%                  functions work.
%
%   The target to work on is determined by on of the following flags:
%
%-     'lib'     - Perform action on the LTFAT C library.
%
%-     'mex'     - Perform action on the mex / oct interfaces.
%
%-     'gpc'     - Perform action on the GPC code for use with MULACLAB
%
%-     'auto'    - Choose automatically which targets to work on based on
%                  the operation system etc. This is the default.
%
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

bp=mfilename('fullpath');
bp=bp(1:end-8);

definput.flags.target={'auto','lib','mex','gpc'};
definput.flags.command={'compile','clean','test'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Remember the current directory.
curdir=pwd;

do_lib  = flags.do_lib || (flags.do_auto && (isoctave || isunix));
do_mex  = flags.do_mex || flags.do_auto;
do_gpc  = flags.do_gpc || (flags.do_auto && ~isoctave);

if isoctave
  extname='oct';
  ext='oct';
else
  extname='mex';
  ext=mexext;
end;

% -------------- Handle cleaning --------------------------------
  
if flags.do_clean
    
  if ~isoctave
    % Delete files permanently (bypass trashcan on Windows/Mac
    % but remember the old state
    oldstate = recycle('off');
  end;
  
  if do_lib
    disp('========= Cleaning libltfat ===============');
    cd([bp,'src']);
    system('make clean');
    disp('Done.');    
  end;
  
  if do_mex
    fprintf('========= Cleaning %s interfaces ==========\n', extname);
    if isoctave
      deletefiles([bp,'oct'],'*.oct');
      deletefiles([bp,'oct'],'*.o');
    else
      deletefiles([bp,'mex'],['*.',mexext]);
    end;
  end;
  
  if do_gpc
    disp('========== Cleaning GPC ==============');
    deletefiles([bp,'thirdparty',filesep,'PolygonClip'],['PolygonClip.',mexext]);
  end;
  
  if ~isoctave
    recycle(oldstate);
  end;  
  
end;

% -------------- Handle compiling  --------------------------------

if flags.do_compile
  if do_lib
    disp('========== Compiling libltfat ==============');
    cd([bp,'src']);
    if isoctave
      [status,output]=system('mkoctfile -p CC');
      system(['make CC=',output]);
    else
      if ispc
        system('make winnomem');
      else
        system('make unixnomem');
      end;
    end;
    disp('Done.');
  end;
  
  if do_mex
    fprintf('========= Compiling %s interfaces ==========\n', extname);
    if compile_ltfat(bp)>1;                
      fprintf('ERROR: The %s interfaces was not built.\n', extname);
    else
      disp('Done.');
    end;
  end;
  
  if do_gpc
    disp('========= Compiling GPC ===============');
    % Compile the PolygonClip interface to GPC for use with mulaclab
    cd([bp,'thirdparty',filesep,'PolygonClip']);
    mex('-I../GPC','PolygonClip.c','../GPC/gpc.c');
    disp('Done.');
  end;
end;

% -------------- Handle testing ---------------------------------------

if flags.do_test
  
  if do_mex
    
    fprintf('========= Testing %s interfaces ==========\n', extname);
    fprintf('1.: Test if comp_pgauss.%s was compiled: ',extname);
    fname=[bp,extname,filesep,'comp_pgauss.',ext];
    if exist(fname,'file')
      disp('SUCCESS.');
    else
      disp('FAILED.');
    end;
    
    fprintf('2.: Test if pgauss executes:              ');
    pgauss(100);
    % If the execution of the script makes it here, we know that pgauss
    % did not crash the system, so we can just print success. Same story
    % with the following entries.
    disp('SUCCESS.');

    fprintf('3.: Test if fftreal executes:             ');
    fftreal(randn(10,1),10);
    disp('SUCCESS.');

    fprintf('4.: Test if dgt executes:                 ');
    dgt(randn(12,1),randn(12,1),3,4);
    disp('SUCCESS.');

    
  end;
  
end;

% Jump back to the original directory.
cd(curdir);


function deletefiles(base,files)

L=dir([base,filesep,files]);
for ii=1:numel(L)
    s=[base,filesep,L(ii).name];
    delete(s);
end;


function status=compile_ltfat(bp)

  uses_lapack = {'comp_gabdual_long','comp_gabtight_long'};
  
% If we exit early, it is because of an error, so set status=1
status=1;

% Determine the name of the ltfat library
s=[bp,'lib',filesep,'libltfat-nomem.a'];

if ispc && ~isoctave
    if strcmp(mexext,'mexw64')
        s=[bp,'mex',filesep,'ltfat.dll'];
    else
        s=[bp,'lib',filesep,'libltfat-nomem.lib'];
    end;
end;

if ~exist(s,'file')
    disp(['The LTFAT C library cannot be found. Please compile ' ...
          'as described in the INSTALL file, or download ' ...
          'a binary release.']);
             
    return;
end;
    
if ispc && ~isoctave
    
    if ~exist([bp,'mex\libfftw3-3.dll'],'file')
        disp(['Please download the FFTW binary release for Windows as ' ...
               'described in the INSTALL file.']);
        return;
    end;
    
    % Create the .lib file for the Matlab FFTW3 lib.
    if ~exist([bp,'mex\libfftw3-3.lib'],'file')
        system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
                bp,'mex\LIBFFTW3-3.DLL"']);
    end;
        
    if ~exist([bp,'mex\libfftw3f-3.dll'],'file')
        disp(['Please download the FFTW binary release for Windows as ' ...
               'described in the INSTALL file.']);
        return;
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
    
    if ~isoctave && strcmp(mexext,'mexa64')
      
      % We don't know how to call LAPACK properly for this platform, so
      % if a mex-file uses LAPACK, skip its compilation
      if any(strcmp(uses_lapack,filename(1:end-endchar-1)))
        disp(['Skipping ',filename]);
        continue
      end;
    end;
    
    if isempty(objdirinfo) || (objdirinfo.datenum<L(ii).datenum)
        
        fprintf('Compiling %s\n',filename);
        
        if isoctave
            % Octave dynamically links to FFTW, Blas and Lapack, so they
            % are not included on the compilation line.
            mkoctfile('-c','oct-memalloc.c');
            mkoctfile('-I../thirdparty',...
                      '-I.','-I../src','-L../src','-L../lib',...
                      filename,'oct-memalloc.o',...
                      '-lltfat-nomem');            
        else
            mex('-c','mex-memalloc.c');
            if ispc
              if strcmp(mexext,'mexw64')
                mex('-I../thirdparty',...
                    '-I.','-I../src','-L../lib','-L../mex',...
                    ['-L',matlabroot,'\extern\lib\win64\microsoft'],...
                    filename,...
                    '-lltfat',...
                    '-lfftw3-3',...
                    '-lfftw3f-3');
                else  
                  mex('-I../thirdparty',...
                      '-I.','-I../src','-L../lib','-L../mex',...
                      ['-L',matlabroot,'\extern\lib\win32\lcc'],...
                      filename,'mex-memalloc.obj',...
                      '-lltfat-nomem',...
                      '-lfftw3-3',...
                      '-lfftw3f-3',...
                      '-lmwlapack',...
                      '-lmwblas');
                end;
                
            else
              
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

status=0;
