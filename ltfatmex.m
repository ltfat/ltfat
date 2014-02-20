function ltfatmex(varargin)
%LTFATMEX   Compile Mex/Oct interfaces
%   Usage:  ltfatmex;
%           ltfatmex(...);
%
%   `ltfatmex` compiles the C backend in order to speed up the execution of
%   the toolbox. The C backend is linked to Matlab and Octave through Mex
%   and Octave C++ interfaces.
%
%   The action of `ltfatmex` is determined by one of the following flags:
%
%     'compile'  Compile stuff. This is the default.
%
%     'clean'    Removes the compiled functions.
%
%     'test'     Run some small tests that verify that the compiled
%                functions work.
%
%   The target to work on is determined by on of the following flags:
%
%     'lib'      Perform action on the LTFAT C library.
%
%     'mex'      Perform action on the mex / oct interfaces.
%
%     'gpc'      Perform action on the GPC code for use with MULACLAB
%
%     'playrec'  Perform action on the playrec code for use with real-time
%                block streaming framework.
%
%     'auto'     Choose automatically which targets to work on based on
%                the operation system etc. This is the default.
%
  
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

bp=mfilename('fullpath');
bp=bp(1:end-length(mfilename));

definput.flags.target={'auto','lib','mex','gpc','playrec','java','blockproc'};
definput.flags.command={'compile','clean','test'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Remember the current directory.
curdir=pwd;

% Compile backend lib?
do_lib  = flags.do_lib || flags.do_auto;
% Compile MEX/OCT interfaces?
do_mex  = flags.do_mex || flags.do_auto;
% Compile MEX PolygonClip.mex... using Generic Polygon Clipper
% (relevant only for the mulaclab, which does not work in Octave)
do_gpc  = flags.do_gpc || (flags.do_auto && ~isoctave);
% Compile MEX playrec.mex... using Portaudio library.
% (relevant only for the bloc processing framework)
do_playrec  = flags.do_playrec || flags.do_blockproc;
% Compile Java classes containing GUI for the bloc proc. framework.
do_java  = flags.do_java || flags.do_blockproc;


if isoctave && flags.do_gpc
    error('%s: Compiling GPC is not relevant for Octave.',upper(mfilename));
end

if isoctave
	extname='oct';
    ext='oct';
else
    extname='mex';
    ext=mexext;
end;

fftw_lib_names = {'fftw3', 'fftw3f' };

% Check if we are on Windows
if ispc
    makefilename='Makefile_mingw';
    make_exe = 'mingw32-make';
    sharedExt = 'dll';
    %fftw_lib_names = {'fftw3', 'fftw3f' };
    % The pre-compiled Octave for Windows comes only in 32bit version (3.6.4)
    % We use different Makefiles
    if isoctave
      makefilename='Makefile_mingwoct';      
    end
end;

% Check if we are on Unix-type system
if isunix
   makefilename='Makefile_unix';
   make_exe = 'make';
   sharedExt = 'so';
end;

% Check if we are on Mac
if ismac
   makefilename='Makefile_mac';
   make_exe = 'make';
   sharedExt = 'dylib';
end;


clear mex;
% -------------- Handle cleaning --------------------------------
if flags.do_clean

  if do_lib
    disp('========= Cleaning libltfat ===============');
    cd([bp,'src']);
    callmake(make_exe,makefilename,'target','clean');
    %[status,result]=system([make_exe, ' -f ',makefilename,' clean']);
    %disp('Done.');    
  end;
  
  if do_mex
    fprintf('========= Cleaning %s interfaces =========\n', upper(extname));
     cd([bp,extname]);
     callmake(make_exe,makefilename,'target','clean','ext',ext);
     %[status,result]=system([make_exe, ' -f ',makefilename,' clean',...
     %                ' EXT=',ext]); 
  end;
  
  if do_gpc
    disp('========= Cleaning GPC ====================');
    cd([bp,'thirdparty',filesep,'PolygonClip']);
    clear mex; 
    callmake(make_exe,makefilename,'target','clean','ext',mexext);
    %[status,result]=system([make_exe, ' -f ',makefilename,' clean',' EXT=',mexext]);
  end;
  
  if do_playrec
     % Use the correct makefile 
     if isoctave
       if ~strcmpi(makefilename(end-2:end),ext)
          makefilename = [makefilename,ext];
       end
    end 
      
    disp('========= Cleaning PLAYREC ================');
    cd([bp,'thirdparty',filesep,'Playrec']);
    clear mex; 
    %[status,result]=system([make_exe, ' -f ',makefilename,' clean',' EXT=',mexext]);
    callmake(make_exe,makefilename,'target','clean','ext',mexext); 
  end;
  
  if do_java
    disp('========= Cleaning JAVA ================');
    cd([bp,'blockproc',filesep,'java']);
    %[status,result]=system([make_exe,' clean']);
    callmake(make_exe,[],'target','clean');
  end;

  cd(curdir);
end;

% -------------- Handle compiling  --------------------------------

if flags.do_compile
  if do_lib
    disp('========= Compiling libltfat ==============');
    cd([bp,'src']);
    clear mex; 
    
    dfftw = ['-l',fftw_lib_names{1}];
    sfftw = ['-l',fftw_lib_names{2}];
    if ispc && ~isoctave
        fftw_lib_found_names = searchfor(bp,fftw_lib_names,sharedExt);
        if ~isempty(fftw_lib_found_names)
           dfftw = ['-l',fftw_lib_found_names{1}(1:end-numel(sharedExt)-1)];
           sfftw = ['-l',fftw_lib_found_names{2}(1:end-numel(sharedExt)-1)];
       end
    end
      % DFFTW and SFFTW are not used in the unix_makefile
      [status,result] = callmake(make_exe,makefilename,'matlabroot','arch',...
                       'dfftw',dfftw,'sfftw',sfftw);
%       [status,result]=system([make_exe, ' -f ',makefilename,...
%                   ' MATLABROOT=','"',matlabroot,'"',...
%                   ' ARCH=',computer('arch'),'"',...
%                   ' DFFTW=',dfftw,'"',...
%                   ' SFFTW=',sfftw]);
      if(~status)
        disp('Done.');
      else
        error('Failed to build LTFAT libs:\n %s',result);
      end
    %end;
  end;
  
  if do_mex
    fprintf('========= Compiling %s interfaces ========\n', upper(extname));
    clear mex; 
    cd([bp,extname]);
    
    dfftw = ['-l',fftw_lib_names{1}];
    sfftw = ['-l',fftw_lib_names{2}];
    if ~isoctave
        fftw_lib_found_names = searchfor(bp,fftw_lib_names,sharedExt);
        if ~isempty(fftw_lib_found_names)
           dfftw = fftw_lib_found_names{1};
           sfftw = fftw_lib_found_names{2};
           if isunix && ~strcmpi(dfftw(1:2),'"/')
              dfftw = ['-l',dfftw(1:end-numel(sharedExt)-1)];
              sfftw = ['-l',sfftw(1:end-numel(sharedExt)-1)];
           end
       end
    end
    
    [status,result] = callmake(make_exe,makefilename,'matlabroot','arch',...
                      'ext',ext,'dfftw',dfftw,'sfftw',sfftw);
%        [status,result]=system([make_exe, ' -f ',makefilename,...
%                             ' MATLABROOT=','"',matlabroot,'"',...
%                             ' EXT=',ext,...
%                             ' ARCH=',computer('arch')]);
    if(~status)
      disp('Done.');
    else
      error('Failed to build %s interfaces: %s \n',upper(extname),result);
    end
  end;
  
  if do_gpc
    disp('========= Compiling GPC ===================');
    % Compile the PolygonClip interface to GPC for use with mulaclab
    cd([bp,'thirdparty',filesep,'PolygonClip']);
    clear mex; 
    [status,result] = callmake(make_exe,makefilename,'matlabroot','arch',...
                      'ext',ext);
%     [status,result]=system([make_exe, ' -f ',makefilename,...
%                      ' MATLABROOT=','"',matlabroot,'"',...
%                      ' EXT=',mexext,...
%                      ' ARCH=',computer('arch')]);
    if(~status)
      disp('Done.');
    else
      error('Failed to build GPC:\n %s',result);
    end
  end;
  if do_playrec
    disp('========= Compiling PLAYREC ===============');
    % Compile the Playrec (interface to portaudio) for the real-time block-
    % stream processing

    if isoctave
        % Because mkoctfile automatically adds -l
       portaudioLib = 'portaudio';
    else
       portaudioLib = '-lportaudio';
    end

       binArchPath = [matlabroot,filesep,'bin',filesep,computer('arch')];
       playrecPath = [bp,'thirdparty',filesep,'Playrec'];
       % Check if portaudio is in thirdparty/Playrec/
       foundPAuser = dir([playrecPath,filesep,'*portaudio*',sharedExt,'*']);
       
       foundPAmatlab = [];
       if ~isoctave
          % Check if portaudio library is present in the Matlab instalation
          foundPAmatlab = dir([binArchPath,filesep,'*portaudio*',sharedExt,'*']);
       end
       
       if ~isempty(foundPAuser)
          if numel(foundPAuser)>1
             error('Ambiguous portaudio libraries in %s. Please leave just one.',playrecPath);
          end
          portaudioLib = foundPAuser(1).name;
          if isunix
              % Need full path
              portaudioLib=[playrecPath,filesep,portaudioLib];
          end
          fprintf('Using %s from thirdparty/Playrec.\n',portaudioLib);
       elseif ~isempty(foundPAmatlab)
          if numel(foundPAmatlab)>1
             if ispc 
                %This should not happen on Windows
                %Use the first one on Linux
                error('Ambiguous portaudio libraries in %s.',binArchPath);
             end
          end
             portaudioLib = foundPAmatlab(1).name;
             if isunix
                portaudioLib = [binArchPath,filesep,portaudioLib];
             end
          
          fprintf('Using %s from Matlab instalation.\n',foundPAmatlab(1).name);
       else
          if ispc && isoctave || ispc
          error(['Portaudio not found. Please download Portaudio http://www.portaudio.com\n',...
                 'and build it as a shared library and copy it to the\n',...
                 '%s directory. \n'],playrecPath);
          elseif ~isoctave
              warning('Portaudio lib should be installed on your system.')
              % Let us hope portaudio is installed on the system
          end
       end


       % Crop library filename
%        if numel(portaudioLib)>3
%          if strcmp(portaudioLib(1:3),'lib')
%            portaudioLib = portaudioLib(4:end);
%          end
%          dotsInPath = strfind(portaudioLib,'.');
%          if ~isempty(dotsInPath) 
%            portaudioLib = portaudioLib(1:dotsInPath(1)-1);
%          end
%        end
    
    if isoctave
       if ~strcmpi(makefilename(end-2:end),ext)
          makefilename = [makefilename,ext];
       end
    end
    
    cd([bp,'thirdparty',filesep,'Playrec']);
    clear mex; 
    [status,result] = callmake(make_exe,makefilename,'matlabroot','arch',...
                      'ext',mexext,'portaudio',portaudioLib);
%     [status,result]=system([make_exe, ' -f ',makefilename,...
%                      ' MATLABROOT=','"',matlabroot,'"',...
%                      ' EXT=',mexext,...
%                      ' PORTAUDIO=',portaudioLib,...
%                      ' ARCH=',computer('arch')]);
    if(~status)
      disp('Done.');
    else
      error('Failed to build PLAYREC:\n %s',result);
    end
  end;
  
  if do_java
    disp('========= Compiling JAVA classes ===================');
    % Compile the JAVA classes
    cd([bp,'blockproc',filesep,'java']);
    clear mex; 
    [status,result] = callmake(make_exe);
    if(~status)
      disp('Done.');
    else
      error('Failed to build JAVA classes:\n %s',result);
    end
  end;
end;

% -------------- Handle testing ---------------------------------------

if flags.do_test
  
  if do_mex
    
    fprintf('========= Testing %s interfaces ==========\n', extname);
    fprintf('1.: Test if comp_pgauss.%s was compiled: ',ext);
    fname=['comp_pgauss.',ext];
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


function status = filesExist(filenames)
   if(~iscell(filenames))
      filenames={filenames};
   end
   for ii=1:length(filenames)
      filename = filenames{ii};
      if(~exist(filename,'file'))
         error('%s: File %s not found.',mfilename,filename);
      end
   end
 
function found_files=searchfor(bp,files,sharedExt)

      found_names = {};
      if ispc 
         for ii=1:numel(files) 
            L = dir([bp,'mex',filesep,'*',files{ii},'*.',sharedExt]);
            if isempty(L)
                error(['%s: %s could not be found in ltfat/mex subdir.',...
                       ' Please download the FFTW dlls and istall them.'],...
                      upper(mfilename),files{ii});
            end
            found_files{ii} = L(1).name;
            fprintf('   ...using %s from ltfat/mex.\n',L(1).name);
         end
      elseif isunix
          
          
          binArchPath = [matlabroot,filesep,'bin',filesep,computer('arch')];
          mexPath = [bp,'mex'];
          
          for ii=1:numel(files)
             % First, try searching in mex dir
             absPath = 1;
             path = mexPath;
             L = dir([mexPath,filesep,'*',files{ii},'*.',sharedExt,'*']); 
             if isempty(L)
                L = dir([binArchPath,filesep,'*',files{ii},'*.',sharedExt,'*']); 
                absPath = 0;
                path = binArchPath;
             end
             
             
             if isempty(L)
                 error('%s: Matlab FFTW libs were not found. Strange.',...
                      upper(mfilename));
             end
             fname = L(1).name;
             if strcmpi(fname(end-2:end),['.',sharedExt]) && ~absPath
             % This is a well behaved case.
             % Just remove lib prefix
             found_files{ii} = fname;
             
             if strcmpi(found_files{ii}(1:3),'lib')
                 found_files{ii} = found_files{ii}(4:end);
             end
             
             
             else    
             % We need a full path here !!
             % Escaped for safety
             found_files{ii} = ['"',path,filesep,fname,'"'];
             end
             fprintf('   ...using %s from Matlab instalation.\n',fname);
          end
          
      end;

   
function [status,result]=callmake(make_exe,makefilename,varargin)
  

  if ~isempty(makefilename) || nargin < 2
     systemCommand = [make_exe, ' -f ',makefilename];
  else
     systemCommand = make_exe; 
  end
  definput.flags.matlabroot={'none','matlabroot'};
  definput.flags.arch={'none','arch'};
  definput.keyvals.ext=[];
  definput.keyvals.dfftw=[];
  definput.keyvals.sfftw=[];
  definput.keyvals.target=[];
  definput.keyvals.portaudio=[];
  [flags,kv]=ltfatarghelper({},definput,varargin);
  
  if flags.do_matlabroot
     systemCommand = [systemCommand, ' MATLABROOT=','"',matlabroot,'"']; 
  end
  
  if flags.do_arch
     systemCommand = [systemCommand, ' ARCH=',computer('arch')]; 
  end
  
  if ~isempty(kv.ext)
     systemCommand = [systemCommand, ' EXT=',kv.ext]; 
  end
  
  if ~isempty(kv.dfftw)
     systemCommand = [systemCommand, ' DFFTW=',kv.dfftw]; 
  end
  
  if ~isempty(kv.sfftw)
     systemCommand = [systemCommand, ' SFFTW=',kv.sfftw]; 
  end
  
  if ~isempty(kv.portaudio)
     systemCommand = [systemCommand, ' PORTAUDIO=',kv.portaudio]; 
  end
  
  if ~isempty(kv.target)
     systemCommand = [systemCommand,' ',kv.target];  
  end

  [status,result]=system(systemCommand);
% 
% function deletefiles(base,files)
% 
% L=dir([base,filesep,files]);
% for ii=1:numel(L)
%     s=[base,filesep,L(ii).name];
%     delete(s);
% end;
% 
% 
% function status=compile_ltfat(bp)
% 
%   uses_lapack = {'comp_gabdual_long','comp_gabtight_long'};
%   
% % If we exit early, it is because of an error, so set status=1
% status=1;
% 
% % Determine the name of the ltfat library
% s=[bp,'lib',filesep,'libltfat-nomem.a'];
% 
% if ispc && ~isoctave
%     if strcmp(mexext,'mexw64')
%         s=[bp,'mex',filesep,'ltfat.dll'];
%     else
%         s=[bp,'lib',filesep,'libltfat-nomem.lib'];
%     end;
% end;
% 
% if ~exist(s,'file')
%     disp(['The LTFAT C library cannot be found. Please compile ' ...
%           'as described in the INSTALL file, or download ' ...
%           'a binary release.']);
%              
%     return;
% end;
%     
% if ispc && ~isoctave
%     
%     if ~exist([bp,'mex\libfftw3-3.dll'],'file')
%         disp(['Please download the FFTW binary release for Windows as ' ...
%                'described in the INSTALL file.']);
%         return;
%     end;
%     
%     % Create the .lib file for the Matlab FFTW3 lib.
%     if ~exist([bp,'mex\libfftw3-3.lib'],'file')
%         system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
%                 bp,'mex\LIBFFTW3-3.DLL"']);
%     end;
%         
%     if ~exist([bp,'mex\libfftw3f-3.dll'],'file')
%         disp(['Please download the FFTW binary release for Windows as ' ...
%                'described in the INSTALL file.']);
%         return;
%     end;
%         
%     if ~exist([bp,'mex\libfftw3f-3.lib'],'file')
%         system(['"',matlabroot,'\sys\lcc\bin\lcc_implib.exe" -u "',...
%                 bp,'mex\libfftw3f-3.dll"']);
%     end;
%         
% end;
%     
% if isoctave
%     cd([bp,'oct']);
%     
%     ext='oct';
%     
%     % Get the list of files.
%     L=dir('*.cc');
%     
%     endchar=2;
% else
%     
%     cd([bp,'mex']);
%     
%     ext=mexext;
%     
%     % Get the list of files.
%     L=dir('comp_*.c');
%     
%     endchar=1;
% end;
% 
% 
% for ii=1:numel(L)
%     filename = L(ii).name;
%     objname  = [filename(1:end-endchar),ext];
%     objdirinfo = dir(objname);
%     
%     % Make-like behaviour: build only the files where the src file is
%     % newer than the object file, or the object file is missing.
%     
%     if ~isoctave && strcmp(mexext,'mexa64')
%       
%       % We don't know how to call LAPACK properly for this platform, so
%       % if a mex-file uses LAPACK, skip its compilation
%       if any(strcmp(uses_lapack,filename(1:end-endchar-1)))
%         disp(['Skipping ',filename]);
%         continue
%       end;
%     end;
%     
%     if isempty(objdirinfo) || (objdirinfo.datenum<L(ii).datenum)
%         
%         fprintf('Compiling %s\n',filename);
%         
%         if isoctave
%             % Octave dynamically links to FFTW, Blas and Lapack, so they
%             % are not included on the compilation line.
%             mkoctfile('-c','oct-memalloc.c');
%             mkoctfile('-I../thirdparty',...
%                       '-I.','-I../src','-L../src','-L../lib',...
%                       filename,'oct-memalloc.o',...
%                       '-lltfat-nomem');            
%         else
%             mex('-c','mex-memalloc.c');
%             if ispc
%               if strcmp(mexext,'mexw64')
%                 mex('-I../thirdparty',...
%                     '-I.','-I../src','-L../lib','-L../mex',...
%                     ['-L',matlabroot,'\extern\lib\win64\microsoft'],...
%                     filename,...
%                     '-lltfat',...
%                     '-lfftw3-3',...
%                     '-lfftw3f-3');
%                 else  
%                   mex('-I../thirdparty',...
%                       '-I.','-I../src','-L../lib','-L../mex',...
%                       ['-L',matlabroot,'\extern\lib\win32\lcc'],...
%                       filename,'mex-memalloc.obj',...
%                       '-lltfat-nomem',...
%                       '-lfftw3-3',...
%                       '-lfftw3f-3',...
%                       '-lmwlapack',...
%                       '-lmwblas');
%                 end;
%                 
%             else
%               
%                 mex('-I../thirdparty',...
%                     '-I.','-I../src','-L../src','-L../lib',...
%                     filename,'mex-memalloc.o',...
%                     '-lltfat-nomem',...
%                     '-lfftw3',...
%                     '-lfftw3f',...
%                     '-llapack',...
%                     '-lblas');
%                 
%             end;
%             
%         end;
%         
%         
%     end;        
%     
% end;
% 
% 
% status=0;
% 
