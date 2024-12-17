function [f,F] = comp_fusyn(F,coeff,varargin)
%FUSYN calculates the fusion frame synthesis operator
%
%   Input parameters:
%     F        fusion frame
%     coeff    coeffcients, given as coeff{1} .. coeff{N}, for each 
%              element in C^L, where L = F.cdim and N= G.Nframes.
%              (can also be a LxN array)
%
%   Output parameters:
%     f   synthesized signal 
%     F   loop back of usion frames
%
%
%   `fusyn` Construct the fusion synthesis operator, i.e. 
%   $\sum \limits_i w_i \pi_{W_i} f_i $, where f_i = coeff{i}.
%   We calculate
%   $ f = \sum \limits_{l = 1}{N} w_l \pi_{W_l} f_l , $, where $(f_l)$ 
%   are the coefficients, $w_l$ the weights, $W_l$ the subspaces
%   (spanned by the local frames).
%
%
%   Optional parameters:
%     'projections'   Turn the projection off (Default: 0), i.e. 
%                     $\sum \limits_i w_i f_i $
%
%     'operator'      can be 'square' or 'frameoper': Square the weights. 
%                     $\sum \limits_i w_i^2 \pi_{W_i} f_i $
%                     We are considering \ell^2(H) as coefficient space, 
%                     so this is equivalent to the frame operator, if we 
%                     use f_i = f.
%
%
% Started: 23.12.2021
% Current: 22.02.2022
%
% Author: P. Balazs
%
% Dependencies: LTFAT http://ltfat.org
%
% Copyright : (c) Acoustics Research Institute, Austrian Academy of
%              Science
%              http://www.kfs.oeaw.ac.at
%
%              Permission is granted to modify and re-distribute this
%              code for non-commercial usage in any manner as long as this
%              notice is preserved. For research please cite the paper [1]. 
%              All standard disclaimers apply.
%
% [1] P. Balazs, M. Shamsabadi, A. Arefijamaal, G. Chardon, "Representation
% of Operators Using Fusion Frames", submitted, 2022


complainif_notvalidframeobj(F,'Fusyn');

if iscell(coeff)
    N_coeff = size(coeff);
    dim_coeff = length(coeff{1});
elseif ismatrix(coeff)
    [dim_coeff,N_coeff] = size(coeff);
    coeff = mat2cell(coeff,dim_coeff,ones(1,N_coeff));
else
    error("The input of fusyn has to be a cell or an array.")
end


definput.flags.projections = {'no_proj', 'proj'};
definput.flags.operator = {'frameoper', 'square'};
[flags,~]=ltfatarghelper({},definput,varargin);
%do_proj = 1;
% Apply the projection, i.e
% $ T_W f_i = sum \limits_i w_i \pi_{W_i} f_i


if ~strcmp(F.type, 'fusion')
    error("Function \'fusyn\' only works for fusion frames");
end
L = size(coeff, 1);
F = comp_checkfudim(F, L);

if isnan(F.cdim)
    error('Local frames have to have the same dimension for the synthesis.');
elseif F.cdim ~= dim_coeff
    error('Dimension of fusion frame and coefficeints do not fit.')
end

if F.Nframes ~= N_coeff
    error("Number of Frames in Fusion Frame and coefficient do not match.")
end


for ii = 1:length(coeff)
     if flags.do_proj
          fff = coeff{ii}; % Only for debugging: Seperation
          ddd = comp_fuana(F,fff,'single',ii);
          if flags.do_square
                ddd = (F.w(ii)).^2*ddd;
          end
          % $= w_i \pi_{W_i} f_i$
     else 
          % $= w_i f_i $
          ddd = F.w(ii).*coeff{ii};
          if flags.do_frameoper
            ddd = {F.w(ii).*ddd};
          end
     end

    if ii ~= 1 
        if length(coeff{ii}) == lll
            f = f{1} + ddd{1};%this is not slower than the cellfun
            if ~iscell(f)
                f={f};
            end
        else 
            error('Dimensions do not match!')
        end
    else
        lll = length(coeff{ii});
        f = ddd;
    end
end