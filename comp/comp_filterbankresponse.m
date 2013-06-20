function gf=comp_filterbankresponse(g,a,L,do_real)

M=numel(g);

if size(a,2)>1
    % G1 is done this way just so that we can determine the data type.
    G1=comp_transferfunction(g{1},L);
    gf=abs(G1).^2*(L/a(1,1)*a(1,2));    
    
    for m=2:M
        gf=gf+abs(comp_transferfunction(g{m},L)).^2*(L/a(m,1)*a(m,2));
    end;
    
else
    % G1 is done this way just so that we can determine the data type.
    G1=comp_transferfunction(g{1},L);
    gf=abs(G1).^2*(L/a(1));    
    
    for m=2:M
        gf=gf+abs(comp_transferfunction(g{m},L)).^2*(L/a(m));
    end;
    
end;
    
if do_real
    gf=gf+involute(gf);   
end;

gf=gf/L;

