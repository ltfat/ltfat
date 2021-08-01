function out = rewrite_struct(in)
    curres = in;
    fn = fieldnames(curres);
    for k=1:numel(fn)
        if( is_function_handle(curres.(fn{k})) )
            curres.(fn{k}) = char(curres.(fn{k}));
        end
    end
    out = curres;
end