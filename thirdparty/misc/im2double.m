% ########################################################################
%
% Copyright (C) 2007-2023 The Octave Project Developers
%
% See the file COPYRIGHT.md in the top-level directory of this
% distribution or <https://octave.org/copyright/>.
%
% This file is part of Octave.
%
% Octave is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <https://www.gnu.org/licenses/>.
%
% ########################################################################
% 
% -*- texinfo -*-
% @deftypefn  {} {@var{dimg} =} im2double (@var{img})
% @deftypefnx {} {@var{dimg} =} im2double (@var{img}, "indexed")
% Convert image to double precision.
%
% The conversion of @var{img} to double precision, is dependent on the type of
% input image.  The following input classes are supported:
%
% @table @samp
% @item uint8, uint16, and int16
% The range of values from the class is scaled to the interval [0 1].
%
% @item logical
% True and false values are assigned a value of 1 and 0 respectively.
%
% @item single
% Values are cast to double.
%
% @item double
% Returns the same image.
%
% @end table
%
% If @var{img} is an indexed image, then the second argument should be the
% string @qcode{"indexed"}.  If so, then @var{img} must either be of floating
% point class, or unsigned integer class and it will simply be cast to double.
% If it is an integer class, an offset of +1 is applied.
%
% @seealso{double}
% @end deftypefn

function dimg = im2double (img, im_type)

if (nargin < 1)
    print_usage ();
end

if (nargin == 1)
    % "normal" (non-indexed) images
    switch (class (img))
    case "uint8",   dimg = double (img) / 255;
    case "uint16",  dimg = double (img) / 65535;
    case "int16",   dimg = (double (img) + 32768) / 65535;
    case "single",  dimg = double (img);
    case "logical", dimg = double (img);
    case "double",  dimg = img;
    otherwise, error ('im2double: IMG is of unsupported class "%s"', class (img));
    end
else
    % indexed images
    if (~ strcmpi (im_type, "indexed"))
      error ('im2double: second input argument must be the string "indexed"');
    elseif (any (isa (img, {"uint8", "uint16"})))
      dimg = double (img) + 1;
    elseif (isfloat (img) || isbool (img))
      dimg = double (img);
    else
      % Technically, it could also be of logical class and we do not
      % enforce positive integers for floating for Matlab compatibility.
      % Still, no need to tell that to the user.
      error (["im2double: if IMG is indexed, then it must be positive " ...
              "integer floating points, or unsigned integer class"]);
    end
end
    
%end
   
  