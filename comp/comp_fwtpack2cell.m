function ccell=comp_fwtpack2cell(F,c)
%COMP_FWTPACK2CELL Change FWT coef. format from pack to cell
%   Usage: ccell=comp_fwtpack2cell(F,c)
%
%   Input parameters:
%         F      : FWT frame object
%         c      : Coefficients in pack format
%   Output parameters:
%         ccell  : Coefficients in cell format
%
%   `comp_fwtpack2cell(F,c)` exctracts individual FWT subbands from 
%   coefficients in packed format *c* as elements of a cell array. *F* must
%   be of type 'fwt' e.g. obtained by `F=frame('fwt',...)` and *c* must
%   be a $Lc \times W$ matrix, obtained by |frana| or |blockana|.
%   
%   The inverse operation is mere c=cell2mat(ccel)
%
%   THE FUNCTION DOES NOT CHECK THE INPUT PARAMETERS IN ANY WAY!
%

w = F.g;
filtNo = numel(w.g);
J = F.J;
subbNo = (filtNo-1)*J+1;
Lc = zeros(subbNo,1);

runPtr = 0;
levelLen = size(c,1)/F.red;

for jj=1:J
     for ff=filtNo:-1:2
        Lc(end-runPtr) = ceil(levelLen/w.a(ff));
        runPtr = runPtr + 1;
     end
     levelLen = ceil(levelLen/w.a(1));
end
Lc(1)=levelLen;

ccell = mat2cell(c,Lc);

% The following does not work for a not being all equal.
%
% filtNo = numel(F.g.g);
% a = F.g.a(1);
% J = F.J;
% 
% subbNo = (filtNo-1)*J+1;
% Lc = zeros(subbNo,1);
% 
% Lc(1:filtNo) = size(c,1)/(1+(filtNo-1)*(a^(J)-1)/(a-1));
% 
% Lc(filtNo+1:end) = kron(Lc(1).*a.^(1:J-1),ones(1,filtNo-1));
% 
% ccell = mat2cell(c,Lc);
