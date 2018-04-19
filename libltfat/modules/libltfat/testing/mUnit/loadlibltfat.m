function loadlibltfat(varargin)
global libltfat_intptrstr;

definput.keyvals.lib='libltfat.so';
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
   
    makecmd = [makecmd, ' MODULE=libltfat'];
    makecmd = [makecmd, ' munit -j12'];
    makecmd = [makecmd, ' MATLABROOT=', matlabroot];
    makecmd = [makecmd, sprintf(' COMPTARGET=%s',flags.comptarget)];
    
    if flags.do_cpp
        makecmd = [makecmd, ' USECPP=1'];
    end
    
    if kv.compiler
        makecmd = [makecmd, sprintf(' CC=%s',kv.compiler)];
    end
    
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
headerpath = [libltfatpath,'build',filesep,'ltfat.h'];
loadlibrary(libpath,headerpath,'mfilename','libltfatprotofile.m');
warning('on');

intbitsize = 8*calllib('libltfat','ltfat_int_size');
libltfat_intptrstr = sprintf('int%dPtr',intbitsize);








