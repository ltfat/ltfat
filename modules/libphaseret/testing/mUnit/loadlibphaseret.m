function loadlibphaseret(varargin)

definput.keyvals.lib='libphaseret.so';
definput.flags.phase={'load','reload'};
[flags,~,lib]=ltfatarghelper({'lib'},definput,varargin);

[~,libname]=fileparts(lib);
    
if libisloaded(libname) 
    if flags.do_reload
        unloadlibrary(libname);
    else
        error('%s: libphaseret is already loaded. Use ''reload'' to force reload.',upper(mfilename));
    end
end

warning('off');
currdir = fileparts(mfilename('fullpath'));
libpath = [currdir, filesep, '..', filesep, '..' , filesep,'build',filesep,lib];
headerpath = [currdir, filesep, '..', filesep, '..' , filesep, 'build',filesep,'phaseret.h'];
loadlibrary(libpath,headerpath,'mfilename','libphaseretprotofile.m');
warning('on');
