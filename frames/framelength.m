function L=framelength(F,Ls);
%FRAMELENGTH  Frame length from signal
%   Usage: L=framelength(F,Ls);
%
%   `framelength(F,Ls)` returns the length of the frame *F*, such that
%   *F* is long enough to expand a signal of length *Ls*.
%
%   If the frame length is longer than the signal length, the signal will be
%   zero-padded by |frana|.
%
%   If instead a set of coefficients are given, call |framelengthcoef|.
%
%   See also: frame, framelengthcoef
  
% Default value, the frame works for all input lengths
L=Ls;

if isfield(F,'framelength')

    L=F.framelength(Ls);

else

    switch(F.type)
      case {'fusion'}
        % This is highly tricky: Get the minimal transform length for each
        % subframe, and set the length as the lcm of that.
        Lsmallest=1;
        for ii=1:F.Nframes
            Lsmallest=lcm(Lsmallest,framelength(F.frames{ii},1));
        end;
        L=ceil(Ls/Lsmallest)*Lsmallest;
        
        % Verify that we did not screw up the assumptions.
        for ii=1:F.Nframes
            if L~=framelength(F.frames{ii},L)
                error(['%s: Cannot determine a frame length. Frame no. %i does ' ...
                       'not support a length of L=%i.'],upper(mfilename),ii,L);
            end;
        end;
      case {'tensor'}
        error(['For the tensor product frame, call framelength for the ' ...
               'individual subframes']);
    end;
    
end;