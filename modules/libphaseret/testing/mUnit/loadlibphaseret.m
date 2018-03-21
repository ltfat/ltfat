function loadlibphaseret(varargin)

definput.keyvals.lib='libphaseret.so';
definput.flags.phase={'load','reload','recompile'};
definput.flags.comptarget={'fulloptim','release','debug'};
definput.flags.verbosity={'quiet','verbose'};
definput.flags.corcpp={'c','cpp'};
definput.keyvals.compiler = [];

[flags,kv,lib]=ltfatarghelper({'lib'},definput,varargin);

[~,libname]=fileparts(lib);
currdir = fileparts(mfilename('fullpath'));
libltfatpath = [currdir, filesep, '..', filesep,'..',filesep,'..',filesep,'..',filesep];
libpath = [libltfatpath, filesep,'build',filesep,lib];
iscompiled = exist(libpath,'file');

makecmd = ['make -C ',libltfatpath];
if ~iscompiled || flags.do_recompile
    [status,result] = system([makecmd, ' clean']);
    if status ~=0, error(result);  end    
   
    makecmd = [makecmd, ' MODULE=%s'];
    makecmd = [makecmd, ' munit -j12'];
    makecmd = [makecmd, ' MATLABROOT=', matlabroot];
    makecmd = [makecmd, sprintf(' COMPTARGET=%s',flags.comptarget)];
    
    if flags.do_cpp
        makecmd = [makecmd, ' USECPP=1'];
    end
    
    if kv.compiler
        makecmd = [makecmd, sprintf(' CC=%s',kv.compiler)];
    end
       
    makecmd_libltfat = sprintf(makecmd,'libltfat');
    if flags.do_verbose
        disp(makecmd_libltfat);
        system(makecmd_libltfat);
    else
        [status,result] = system(makecmd_libltfat);
        if status ~=0, error(result);  end
    end
        
    makecmd_libphaseret = sprintf(makecmd,'libphaseret');
    if flags.do_verbose
        disp(makecmd_libphaseret);
        system(makecmd_libphaseret);
    else
        [status,result] = system(makecmd_libphaseret);
        if status ~=0, error(result);  end
    end
end
    
if libisloaded(libname) 
    if flags.do_reload || flags.do_recompile 
        unloadlibrary(libname);
    else
        error('%s: libphaseret is already loaded. Use ''reload'' to force reload.',upper(mfilename));
    end
end

warning('off');
headerpath = [libltfatpath,'build',filesep,'phaseret.h'];
headerpath_libltfat = [libltfatpath,'build',filesep,'ltfat.h'];
loadlibrary(libpath,headerpath,'mfilename','libphaseretprotofile.m','addheader',headerpath_libltfat);
warning('on');
