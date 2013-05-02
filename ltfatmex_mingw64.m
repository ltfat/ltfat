function buildltfat(varargin)
% 1) Download and unpack 64bit binaries of the FFTW library from http://www.fftw.org/install/windows.html. Copy all *.dll files to the ltfat/mex directory.
% 2) Download and install 64bit TDM-GCC compiler suite from http://tdm-gcc.tdragon.net/download. Ensure that the [MinGW64]/bin directory is in the PATH, where [MinGW64] stands for instalation direcory.
% 3) Run this script
%% 
fftw_lib_names = {'libfftw3-3', 'libfftw3f-3' };
make_exe = 'mingw32-make';
makefilename='Makefile_mingw64';
%%

definput.flags.command={'compile','clean'};
[flags,kv]=ltfatarghelper({},definput,varargin);



bp=mfilename('fullpath');
bp=bp(1:end-length(mfilename));
curdir=pwd;

if flags.do_clean
    clear functions;
    oldstate = recycle('off');
    cd([bp,'thirdparty',filesep,'PolygonClip']);
    [status,result]=system([make_exe, ' -f ',makefilename,' clean',' MEXEXT=',mexext]);
    cd([bp,'src']);
    [status,result]=system([make_exe, ' -f ',makefilename,' clean']);
    cd([bp,'mex']);
    [status,result]=system([make_exe, ' -f ',makefilename,' clean',...
                    ' MEXEXT=',mexext]);
%     for ii=1:length(fftw_lib_names)
%        filename = [fftw_lib_names{ii},'.lib'];
%        delete([bp,'mex\',filename]);
%     end
          
    cd(curdir);
    return;
end

clear functions;
% Compile PolygonClip
disp('========= Compiling GPC ====================================');
cd([bp,'thirdparty',filesep,'PolygonClip']);
[status,result]=system([make_exe, ' -f ',makefilename,...
                  ' MATLABROOT=','"',matlabroot,'"',...
                  ' MEXEXT=',mexext]);
if(~status)
   disp('Done.');
else
   error('Failed to build GPC:\n %s',result);
end

disp('========= Checking for FFTW ================================');
cd(bp);
% check ltfat/mex for DLLs and DEFs
for ii=1:length(fftw_lib_names)
%     filename = [fftw_lib_names{ii},'.def'];
%     if(~exist([bp,'mex\',filename],'file'))
%          error('%s: Definition file %s not found.',mfilename,filename);
%     end
    
    filename = [fftw_lib_names{ii},'.dll'];
    if(~exist([bp,'mex\',filename],'file'))
         error('%s: DLL file %s not found.',mfilename,filename);
    end
end

% Create import .lib files for the Matlab FFTW3 lib. Is probably not
% necessary, as gcc can link directly against dlls.
% for ii=1:length(fftw_lib_names)
%     filename = fftw_lib_names{ii};
%     if ~exist([bp,'mex/',filename,'.lib'],'file')
%         if(system(sprintf('dlltool -d %s -l %s',['mex/',filename,'.def'],['mex/',filename,'.lib'])))
%             error('FAILED. Import library %s.lib was not created.',filename);
%         end
%     end;
% end
 disp('Done.');

disp('========= Building LTFAT ===================================');
% build ltfat.dll
cd([bp,'src']);
[status,result]=system([make_exe, ' -f ',makefilename,...
                  ' MATLABROOT=','"',matlabroot,'"']);
if(~status)
   disp('Done.');
else
   error('Failed to build LTFAT dlls:\n %s',result);
end

disp('========= Building MEX interfaces ===========================');
cd([bp,'mex']);
[status,result]=system([make_exe, ' -f ',makefilename,...
                  ' MATLABROOT=','"',matlabroot,'"',...
                  ' MEXEXT=',mexext]);
if(~status)
   disp('Done.');
else
   error('Failed to build MEX interfaces: %s \n',result);
end

cd(curdir);





