function ltfatmex(varargin)
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

defnopos.flags.target={'auto','lib','mex','oct','gpc'};
defnopos.flags.command={'compile','clean'};
[flags,kv]=ltfatarghelper({},defnopos,varargin);

% Remember the current directory.
curdir=pwd;

do_lib  = flags.do_lib || (flags.do_auto && isoctave);
do_oct  = flags.do_oct || (flags.do_auto && isoctave);
do_mex  = flags.do_mex || (flags.do_auto && ~isoctave);
do_gpc  = flags.do_gpc || (flags.do_auto && ~isoctave);

if flags.do_clean
    
    if ~isoctave
        % Delete files permanently (bypass trashcan on Windows/Mac
        % but remember the old state
        oldstate = recycle('off');
    end;
        
    if do_oct
        deletefiles([bp,'oct'],'*.oct');
        deletefiles([bp,'oct'],'*.o');
    end;
    
    if do_mex
        deletefiles([bp,'mex'],['*.',mexext]);
    end;
    
    if do_gpc
        deletefiles([bp,'thirdparty',filesep,'PolygonClip'],
        ['PolygonClip.',mexext]);
    end;
    
    if ~isoctave
        recycle(oldstate);
    end;
    
    return;
    
end;

if do_lib
    s=['cd ',bp,'src; make'];
    s
    system(s);    
end;

if do_oct
    if compile_ltfat(bp)>1;
        
        if isoctave
            extname='oct';
        else
            extname='mex';
        end;
        
        fprintf('The %s interfaces was not built.', extname);
    end;
end;

    

if do_gpc &&  ~isoctave
    if exist([bp,filesep,'thirdparty',filesep,'GPC',filesep,'gpc.c'], ...
             'file')
        % Compile the PolygonClip interface to GPC for use with mulaclab
        cd([bp,filesep,'thirdparty',filesep,'PolygonClip']);
        mex('-I../GPC','PolygonClip.c','../GPC/gpc.c');
    else
        error('Unable to compile GPC');
        
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

% If we exit early, it is because of an error, so set status=1
status=1;

% Determine the name of the ltfat library
s=[bp,'lib',filesep,'libltfat-nomem.a'];

if ispc || ~isoctave    
    s=[bp,'lib',filesep,'libltfat-nomem.lib'];
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

status=0;