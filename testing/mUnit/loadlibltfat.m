function loadlibltfat(varargin)

definput.keyvals.lib='libltfat.so';
definput.flags.phase={'load','reload','recompile'};
definput.flags.comptarget={'fulloptim','release','debug'};
definput.flags.verbosity={'quiet','verbose'};
[flags,~,lib]=ltfatarghelper({'lib'},definput,varargin);

[~,libname]=fileparts(lib);
currdir = fileparts(mfilename('fullpath'));
libltfatpath = [currdir, filesep, '..', filesep,'..'];
libpath = [libltfatpath, filesep,'build',filesep,lib];
iscompiled = exist(libpath,'file');

makecmd = ['make -C ',libltfatpath];
if ~iscompiled || flags.do_recompile
    makecmd = [makecmd, ' munit -j12'];
    makecmd = [makecmd, ' MATLABROOT=', matlabroot];
    makecmd = [makecmd, sprintf(' COMPTARGET=%s',flags.comptarget)];
    
    if flags.do_verbose
        disp(makecmd);
        system(makecmd);
    else
        [status,result] = system(makecmd);
        if status ~=0, error(result);  end
    end
end
    
if libisloaded(libname) 
    if flags.do_reload || flags.do_recompile 
        unloadlibrary(libname);
    else
        error('%s: libltfat is already loaded. Use ''reload'' to force reload.',upper(mfilename));
    end
end

warning('off');
headerpath = [currdir, filesep,'..', filesep, '..', filesep,'build',filesep,'ltfat.h'];
loadlibrary(libpath,headerpath,'mfilename','libltfatprotofile.m');
warning('on');






