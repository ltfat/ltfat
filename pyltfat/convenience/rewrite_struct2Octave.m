function [out, ret] = rewrite_struct2Octave(in)
    curres = in;
    fn = fieldnames(curres);
    for k=1:numel(fn)
        if isfield(curres, fn(k)) && ischar(curres.(fn{k})) && strcmp(curres.(fn{k})(1), '@')
            curres.(fn{k}) = str2func(curres.(fn{k}));
        end
    end
    out = curres;
    ret = 1;
end